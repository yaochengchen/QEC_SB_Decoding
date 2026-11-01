from dataclasses import dataclass, asdict
import json, time
from pathlib import Path
import numpy as np
import pandas as pd
from .rotated_code import build_rotated_surface_code, sample_depolarizing
from .oracle import build_d3_oracle, syndrome_key, independent_decode
from .qubo import build_depolarizing_qubo, qubo_to_zero_field_ising, warm_binary_from_correction
from .sb import SBParams, run_sb, binary_from_augmented_spins

@dataclass
class BenchConfig:
    distance: int = 3
    p: float = 0.08
    shots: int = 240
    seed: int = 20261008
    penalty_mult: float = 5.0
    replicas: int = 64
    steps: int = 800
    dt: float = 0.5
    gb_control_A: float = 0.2


def _eval_candidates(vv, code, sx, sz, S):
    n = code.n
    ex = vv[:, :n]
    ez = vv[:, n:2*n]
    feas = np.all(((ex @ code.Hz.T) & 1) == sx[None, :], axis=1)
    feas &= np.all(((ez @ code.Hx.T) & 1) == sz[None, :], axis=1)
    if not np.any(feas):
        return None
    ids = np.flatnonzero(feas)
    wt = np.count_nonzero(ex[ids] | ez[ids], axis=1)
    qenergy = np.einsum('bi,ij,bj->b', vv[ids], S, vv[ids])
    # QUBO objective equals Pauli weight (up to w) on feasible states.
    order = np.lexsort((qenergy, wt))
    k = ids[order[0]]
    lx = int((ez[k] @ code.logical_x_support) & 1)
    lz = int((ex[k] @ code.logical_z_support) & 1)
    return int(wt[order[0]]), np.asarray([lx, lz], dtype=np.uint8), int(k)


def wilson(k, n, z=1.96):
    if n == 0:
        return (float('nan'), float('nan'))
    p = k / n
    den = 1 + z*z/n
    ctr = (p + z*z/(2*n))/den
    half = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n))/den
    return ctr-half, ctr+half


def run_benchmark(outdir, cfg=BenchConfig()):
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    code = build_rotated_surface_code(cfg.distance)
    if cfg.distance != 3:
        raise ValueError("Current exact ground-state benchmark is intentionally d=3.")
    t0 = time.time()
    oracle = build_d3_oracle(code)
    samples = sample_depolarizing(code, cfg.p, cfg.shots, cfg.seed)

    method_names = ["bSB", "dSB", "GbSB", "MWPM->bSB", "MWPM->dSB", "MWPM->GbSB"]
    stats = {m: {"feasible":0,"ground":0,"logical":0,"logical_not_ground":0,
                 "ground_but_logical_wrong":0,"improves_incumbent_energy":0,
                 "runtime_s":0.0} for m in method_names}
    baseline = {"independent_matching_logical":0, "map_config_logical":0, "exact_mld_logical":0,
                "matching_is_map_ground":0}
    rows = []

    params = SBParams(replicas=cfg.replicas, steps=cfg.steps, dt=cfg.dt,
                      warm_strength=0.08, gb_control_A=cfg.gb_control_A)

    for sh in range(cfg.shots):
        sx, sz = samples["sx"][sh], samples["sz"][sh]
        truth = samples["logical"][sh]
        ex0, ez0, match_logical = independent_decode(code, oracle, sx, sz)
        base_wt = int(np.count_nonzero(ex0 | ez0))
        baseline["independent_matching_logical"] += int(np.array_equal(match_logical, truth))

        key = syndrome_key(sx, sz)
        opt_wt = oracle.min_pauli_weight[key]
        map_ex, map_ez = oracle.map_config[key]
        map_logical = np.asarray([
            int((map_ez @ code.logical_x_support) & 1),
            int((map_ex @ code.logical_z_support) & 1)], dtype=np.uint8)
        baseline["map_config_logical"] += int(np.array_equal(map_logical, truth))
        baseline["matching_is_map_ground"] += int(base_wt == opt_wt)
        mld_logical, _ = oracle.exact_mld_class(key, cfg.p)
        baseline["exact_mld_logical"] += int(np.array_equal(mld_logical, truth))

        model = build_depolarizing_qubo(code, sx, sz, cfg.p, cfg.penalty_mult)
        J, _ = qubo_to_zero_field_ising(model.S)
        warm_v = warm_binary_from_correction(code, model, ex0, ez0, sx, sz)
        warm_spin = np.concatenate([1 - 2 * warm_v.astype(float), [1.0]])

        for mi, method in enumerate(method_names):
            variant = method.split("->")[-1].lower()
            warm = method.startswith("MWPM")
            ts = time.perf_counter()
            spins = run_sb(J, variant, params, seed=cfg.seed + sh*37 + mi*100003,
                           warm_spin=warm_spin if warm else None)
            vv = binary_from_augmented_spins(spins)
            # Keep the MWPM incumbent in the warm-start pool: refinement must not
            # discard the valid seed merely because the nonlinear trajectory moved.
            if warm:
                vv = np.vstack([warm_v[None, :], vv])
            ans = _eval_candidates(vv, code, sx, sz, model.S)
            stats[method]["runtime_s"] += time.perf_counter() - ts
            if ans is None:
                rows.append({"shot":sh,"method":method,"feasible":0})
                continue
            wt, logical, _ = ans
            logical_ok = bool(np.array_equal(logical, truth))
            ground = wt == opt_wt
            st = stats[method]
            st["feasible"] += 1
            st["ground"] += int(ground)
            st["logical"] += int(logical_ok)
            st["logical_not_ground"] += int(logical_ok and not ground)
            st["ground_but_logical_wrong"] += int(ground and not logical_ok)
            st["improves_incumbent_energy"] += int(wt < base_wt)
            rows.append({"shot":sh,"method":method,"feasible":1,"weight":wt,
                         "opt_weight":opt_wt,"ground":int(ground),"logical_ok":int(logical_ok),
                         "logical_not_ground":int(logical_ok and not ground),
                         "ground_but_logical_wrong":int(ground and not logical_ok),
                         "incumbent_weight":base_wt})

    detail = pd.DataFrame(rows)
    detail.to_csv(outdir / "shot_level.csv", index=False)
    summary_rows = []
    for method, st in stats.items():
        n = cfg.shots
        lo, hi = wilson(st["logical"], n)
        summary_rows.append({
            "method": method,
            "feasible_rate": st["feasible"]/n,
            "ground_state_success": st["ground"]/n,
            "logical_success": st["logical"]/n,
            "logical_success_ci95_lo": lo,
            "logical_success_ci95_hi": hi,
            "logical_not_ground_rate": st["logical_not_ground"]/n,
            "ground_but_logical_wrong_rate": st["ground_but_logical_wrong"]/n,
            "energy_improvement_over_matching_rate": st["improves_incumbent_energy"]/n,
            "mean_solver_time_ms_per_shot": 1000*st["runtime_s"]/n,
        })
    base_rows = []
    for name, count in baseline.items():
        lo, hi = wilson(count, cfg.shots)
        base_rows.append({"baseline":name,"success_rate":count/cfg.shots,"ci95_lo":lo,"ci95_hi":hi})
    summary = pd.DataFrame(summary_rows)
    bases = pd.DataFrame(base_rows)
    summary.to_csv(outdir / "summary.csv", index=False)
    bases.to_csv(outdir / "baselines.csv", index=False)
    meta = asdict(cfg)
    meta.update({"elapsed_s":time.time()-t0,"max_exact_map_weight_over_all_syndromes":max(oracle.min_pauli_weight.values())})
    (outdir / "metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    return summary, bases, detail, meta
