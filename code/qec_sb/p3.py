from dataclasses import dataclass, asdict
from itertools import combinations
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from .rotated_code import build_rotated_surface_code, sample_depolarizing
from .logical_classes import logical_class as _logical_class, logical_class_seed as _class_seed
from .qubo import build_depolarizing_qubo, qubo_to_zero_field_ising, warm_binary_from_correction
from .sb import SBParams, run_sb, binary_from_augmented_spins


@dataclass
class P3Config:
    distance: int = 3
    p: float = 0.08
    shots: int = 10000
    seed: int = 20261010
    penalty_mult: float = 5.0
    replicas: int = 64
    steps: int = 600
    dts: tuple = (0.35, 0.50, 0.75)
    top_m: int = 1
    orbit_radius: int = 1
    warm_strength: float = 0.08


def _stable_argmax(scores, tol=1e-15):
    """Return the lowest-index maximum, using a scale-relative tie tolerance."""
    scores = np.asarray(scores, dtype=float)
    mx = float(np.max(scores))
    return int(np.flatnonzero(np.isclose(scores, mx, rtol=tol, atol=0.0))[0])


def _stabilizer_masks(code, radius):
    """Return physical (X|Z) masks for all stabilizer products up to Hamming radius."""
    gens = []
    z0 = np.zeros(code.n, dtype=np.uint8)
    for row in code.Hx:
        gens.append(np.concatenate([row.astype(np.uint8), z0]))
    for row in code.Hz:
        gens.append(np.concatenate([z0, row.astype(np.uint8)]))
    gens = np.asarray(gens, dtype=np.uint8)

    masks = [np.zeros(2 * code.n, dtype=np.uint8)]
    for r in range(1, int(radius) + 1):
        for comb in combinations(range(len(gens)), r):
            mask = np.zeros(2 * code.n, dtype=np.uint8)
            for j in comb:
                mask ^= gens[j]
            masks.append(mask)
    return np.asarray(masks, dtype=np.uint8)


def collect_dsb_candidates(code, oracle, sx, sz, p, seed, cfg: P3Config):
    """Harvest a diverse syndrome-feasible dSB pool.

    For each dt we mix matching-warm and random-start replicas. Only distinct
    physical Pauli configurations are retained; QUBO slack bits are discarded.
    """
    from .oracle import independent_decode

    ex0, ez0, _ = independent_decode(code, oracle, sx, sz)
    model = build_depolarizing_qubo(code, sx, sz, p, cfg.penalty_mult)
    J, _ = qubo_to_zero_field_ising(model.S)
    warm_v = warm_binary_from_correction(code, model, ex0, ez0, sx, sz)
    warm_spin = np.concatenate([1 - 2 * warm_v.astype(float), [1.0]])

    pools = [warm_v[None, :]]
    for j, dt in enumerate(cfg.dts):
        for warm in (True, False):
            params = SBParams(
                replicas=cfg.replicas,
                steps=cfg.steps,
                dt=float(dt),
                warm_strength=cfg.warm_strength,
            )
            spins = run_sb(
                J,
                "dsb",
                params,
                seed=int(seed + 97 * j + (10000 if warm else 0)),
                warm_spin=warm_spin if warm else None,
            )
            pools.append(binary_from_augmented_spins(spins))

    vv = np.vstack(pools)
    n = code.n
    ex = vv[:, :n]
    ez = vv[:, n : 2 * n]
    feasible = np.all(((ex @ code.Hz.T) & 1) == sx[None, :], axis=1)
    feasible &= np.all(((ez @ code.Hx.T) & 1) == sz[None, :], axis=1)
    phys = np.concatenate([ex[feasible], ez[feasible]], axis=1)
    if len(phys) == 0:
        phys = np.concatenate([ex0, ez0])[None, :]
    phys = np.unique(phys.astype(np.uint8), axis=0)
    return phys, ex0, ez0, int(feasible.sum()), int(len(vv))


def naive_truncated_class(pool, code, p):
    """Naive P3 baseline: sum physical Boltzmann weights of unique raw SB samples."""
    n = code.n
    ex, ez = pool[:, :n], pool[:, n:]
    wt = np.count_nonzero(ex | ez, axis=1)
    cls = ((ez @ code.logical_x_support) & 1) * 2 + ((ex @ code.logical_z_support) & 1)
    ratio = (p / 3.0) / (1.0 - p)
    scores = np.zeros(4, dtype=float)
    for c in range(4):
        scores[c] = float(np.sum(ratio ** wt[cls == c]))
    return _stable_argmax(scores), scores


def sb_stabilizer_orbit_class(pool, code, ex0, ez0, p, top_m=1, orbit_radius=1):
    """P3 QEC-aware class scoring via SB seeds + local stabilizer-orbit expansion.

    Per logical class:
      1) keep top-M lowest-Pauli-weight SB seeds,
      2) add a guaranteed class-feasible matching-derived seed,
      3) multiply each seed by stabilizer products up to `orbit_radius`,
      4) deduplicate and sum exact depolarizing relative probabilities.

    The local orbit expansion explicitly restores nearby degenerate errors that a
    biased optimizer sampler can miss. radius=1 costs only 1 + (# stabilizers)
    neighbors per seed and is scalable, unlike full stabilizer enumeration.
    """
    n = code.n
    ex, ez = pool[:, :n], pool[:, n:]
    wt = np.count_nonzero(ex | ez, axis=1)
    cls = ((ez @ code.logical_x_support) & 1) * 2 + ((ex @ code.logical_z_support) & 1)
    masks = _stabilizer_masks(code, orbit_radius)
    ratio = (p / 3.0) / (1.0 - p)

    scores = np.zeros(4, dtype=float)
    expanded_counts = np.zeros(4, dtype=int)
    for c in range(4):
        ids = np.flatnonzero(cls == c)
        seeds = []
        if len(ids):
            order = ids[np.argsort(wt[ids], kind="stable")[: int(top_m)]]
            seeds.extend(pool[order])
        cex, cez = _class_seed(code, ex0, ez0, c)
        seeds.append(np.concatenate([cex, cez]))

        expanded = np.vstack([seed[None, :] ^ masks for seed in seeds])
        expanded = np.unique(expanded.astype(np.uint8), axis=0)
        eex, eez = expanded[:, :n], expanded[:, n:]
        ewt = np.count_nonzero(eex | eez, axis=1)
        scores[c] = float(np.sum(ratio ** ewt))
        expanded_counts[c] = len(expanded)

    return _stable_argmax(scores), scores, expanded_counts


def _paired_diff_ci(a, b, z=1.96):
    a = np.asarray(a, dtype=np.int8)
    b = np.asarray(b, dtype=np.int8)
    d = a - b
    diff = float(d.mean())
    se = float(d.std(ddof=1) / np.sqrt(len(d))) if len(d) > 1 else float("nan")
    return diff, diff - z * se, diff + z * se, int(np.sum((a == 1) & (b == 0))), int(np.sum((a == 0) & (b == 1)))


def run_p3_confirmation(outdir, cfg=P3Config()):
    from .oracle import build_d3_oracle, syndrome_key, independent_decode
    from .benchmark import wilson

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if cfg.distance != 3:
        raise ValueError("P3 exact-oracle confirmation is intentionally d=3.")

    t0 = time.time()
    code = build_rotated_surface_code(cfg.distance)
    oracle = build_d3_oracle(code)
    samples = sample_depolarizing(code, cfg.p, cfg.shots, cfg.seed)
    keys = [syndrome_key(samples["sx"][i], samples["sz"][i]) for i in range(cfg.shots)]

    representative = {}
    for i, key in enumerate(keys):
        representative.setdefault(key, i)

    pools = {}
    pool_rows = []
    for key, i in representative.items():
        pool, ex0, ez0, feasible_count, raw_count = collect_dsb_candidates(
            code,
            oracle,
            samples["sx"][i],
            samples["sz"][i],
            cfg.p,
            cfg.seed + 1000 + int(key),
            cfg,
        )
        pools[key] = (pool, ex0, ez0)
        pool_rows.append(
            {
                "syndrome_key": key,
                "raw_trajectories_plus_incumbent": raw_count,
                "feasible_raw": feasible_count,
                "unique_feasible_physical": len(pool),
            }
        )

    methods = ["matching", "single_MAP_config", "naive_SB_partition", "P3_SB_plus_SOE", "exact_logical_MLD"]
    correct = {m: np.zeros(cfg.shots, dtype=np.uint8) for m in methods}
    p3_agrees_oracle = np.zeros(cfg.shots, dtype=np.uint8)
    naive_agrees_oracle = np.zeros(cfg.shots, dtype=np.uint8)
    shot_rows = []

    for i, key in enumerate(keys):
        sx, sz = samples["sx"][i], samples["sz"][i]
        truth = int(samples["logical"][i, 0]) * 2 + int(samples["logical"][i, 1])
        pool, ex0, ez0 = pools[key]

        _, _, match_log = independent_decode(code, oracle, sx, sz)
        match_cls = int(match_log[0]) * 2 + int(match_log[1])

        map_ex, map_ez = oracle.map_config[key]
        map_cls = _logical_class(code, map_ex, map_ez)

        mld_log, _ = oracle.exact_mld_class(key, cfg.p)
        mld_cls = int(mld_log[0]) * 2 + int(mld_log[1])

        naive_cls, _ = naive_truncated_class(pool, code, cfg.p)
        p3_cls, _, expanded_counts = sb_stabilizer_orbit_class(
            pool,
            code,
            ex0,
            ez0,
            cfg.p,
            top_m=cfg.top_m,
            orbit_radius=cfg.orbit_radius,
        )

        preds = {
            "matching": match_cls,
            "single_MAP_config": map_cls,
            "naive_SB_partition": naive_cls,
            "P3_SB_plus_SOE": p3_cls,
            "exact_logical_MLD": mld_cls,
        }
        for name, pred in preds.items():
            correct[name][i] = int(pred == truth)
        p3_agrees_oracle[i] = int(p3_cls == mld_cls)
        naive_agrees_oracle[i] = int(naive_cls == mld_cls)
        shot_rows.append(
            {
                "shot": i,
                "syndrome_key": key,
                "truth_class": truth,
                **{f"pred_{k}": v for k, v in preds.items()},
                **{f"correct_{k}": int(v == truth) for k, v in preds.items()},
                "p3_expanded_candidates_total": int(expanded_counts.sum()),
            }
        )

    summary_rows = []
    for name in methods:
        k = int(correct[name].sum())
        lo, hi = wilson(k, cfg.shots)
        summary_rows.append({"method": name, "logical_success": k / cfg.shots, "ci95_lo": lo, "ci95_hi": hi})
    summary = pd.DataFrame(summary_rows)

    comparisons = []
    for a, b in [
        ("P3_SB_plus_SOE", "matching"),
        ("P3_SB_plus_SOE", "single_MAP_config"),
        ("P3_SB_plus_SOE", "naive_SB_partition"),
        ("exact_logical_MLD", "P3_SB_plus_SOE"),
    ]:
        diff, lo, hi, a_only, b_only = _paired_diff_ci(correct[a], correct[b])
        comparisons.append(
            {
                "method_a": a,
                "method_b": b,
                "accuracy_diff_a_minus_b": diff,
                "paired_ci95_lo": lo,
                "paired_ci95_hi": hi,
                "a_correct_b_wrong": a_only,
                "a_wrong_b_correct": b_only,
            }
        )

    pd.DataFrame(shot_rows).to_csv(outdir / "shot_level.csv", index=False)
    summary.to_csv(outdir / "summary.csv", index=False)
    pd.DataFrame(comparisons).to_csv(outdir / "paired_comparisons.csv", index=False)
    pool_df = pd.DataFrame(pool_rows)
    pool_df.to_csv(outdir / "candidate_pool_stats.csv", index=False)

    meta = asdict(cfg)
    meta.update(
        {
            "unique_syndromes": len(representative),
            "mean_unique_feasible_candidates_per_syndrome": float(pool_df["unique_feasible_physical"].mean()),
            "p3_exact_class_agreement": float(p3_agrees_oracle.mean()),
            "naive_exact_class_agreement": float(naive_agrees_oracle.mean()),
            "elapsed_s": time.time() - t0,
        }
    )
    (outdir / "metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    return summary, pd.DataFrame(comparisons), pool_df, meta
