from dataclasses import dataclass, asdict
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from .rotated_code import build_rotated_surface_code, sample_depolarizing
from .p3 import _class_seed, _logical_class, _paired_diff_ci, _stable_argmax
from .cp_sb import (
    build_stabilizer_hamiltonian,
    run_constraint_preserving_sb,
    CPSBParams,
    stabilizer_hamming_neighborhood,
)


@dataclass
class CPAdaptiveConfig:
    distance: int = 3
    p: float = 0.08
    shots: int = 10000
    seed: int = 20261015
    solver_seed: int = 20261014
    replicas: int = 8
    steps: int = 120
    dt: float = 0.45
    coupling_scale: float = 0.70
    warm_strength: float = 0.08
    # Candidate diversity: warm + random trajectories for each scale multiplier.
    coupling_multipliers: tuple = (1.0,)
    # Adaptive SOE stages (top_m, Hamming radius in stabilizer variables).
    stages: tuple = ((1, 1), (1, 8))
    # Disable confidence-only stopping by default; use a 5% top/runner-up gap.
    confidence_threshold: float = 1.10
    min_ratio_threshold: float = 1.05


def syndrome_from_key(key: int, n_sx: int, n_sz: int):
    n = n_sx + n_sz
    bits = np.asarray([(int(key) >> (n - 1 - i)) & 1 for i in range(n)], dtype=np.uint8)
    return bits[:n_sx], bits[n_sx:]


def _class_masses_exact(oracle, key, p):
    counts = oracle.class_weight_counts[key]
    r = (p / 3.0) / (1.0 - p)
    return counts @ (r ** np.arange(counts.shape[1]))


def collect_cp_class_pools(code, oracle, sx, sz, base_seed: int, cfg: CPAdaptiveConfig):
    """Collect stabilizer-variable CP-SB candidates separately in all 4 logical classes.

    Every candidate is feasible by construction. Returned pools are generator-bit
    coordinates, making local stabilizer-orbit expansion a Hamming neighborhood.
    """
    from .oracle import independent_decode

    ex0, ez0, _ = independent_decode(code, oracle, sx, sz)
    pools = {}
    hams = {}
    stats = {}
    for c in range(4):
        cex, cez = _class_seed(code, ex0, ez0, c)
        ham = build_stabilizer_hamiltonian(code, cex, cez)
        hams[c] = ham
        raw = [np.zeros((1, ham.n_spins), dtype=np.uint8)]
        for j, mult in enumerate(cfg.coupling_multipliers):
            for warm in (True, False):
                par = CPSBParams(
                    replicas=cfg.replicas,
                    steps=cfg.steps,
                    dt=cfg.dt,
                    coupling_scale=cfg.coupling_scale * float(mult),
                    warm_strength=cfg.warm_strength,
                )
                spins = run_constraint_preserving_sb(
                    ham,
                    par,
                    seed=int(base_seed + c * 100003 + j * 997 + (500000 if warm else 0)),
                    warm_bits=np.zeros(ham.n_spins, dtype=np.uint8) if warm else None,
                )
                raw.append((spins < 0).astype(np.uint8))
        u = np.unique(np.vstack(raw), axis=0)
        pools[c] = u
        phys = ham.physical_from_bits(u)
        wt = np.count_nonzero(phys[:, :code.n] | phys[:, code.n:], axis=1)
        stats[c] = {
            "unique_cp_candidates": int(len(u)),
            "best_weight": int(wt.min()),
            "mean_weight": float(wt.mean()),
        }
    return pools, hams, ex0, ez0, stats


def _score_one_class(pool_u, ham, p, top_m, radius):
    phys = ham.physical_from_bits(pool_u)
    n = ham.n_data
    wt = np.count_nonzero(phys[:, :n] | phys[:, n:], axis=1)
    order = np.argsort(wt, kind="stable")[: max(1, int(top_m))]
    expanded = []
    for idx in order:
        expanded.append(stabilizer_hamming_neighborhood(pool_u[idx], ham.n_spins, radius))
    uu = np.unique(np.vstack(expanded), axis=0)
    e = ham.physical_from_bits(uu)
    ewt = np.count_nonzero(e[:, :n] | e[:, n:], axis=1)
    ratio = (p / 3.0) / (1.0 - p)
    score = float(np.sum(ratio ** ewt))
    return score, int(len(uu)), int(ewt.min()), uu


def fixed_soe_prediction(pools, hams, p, top_m=1, radius=1):
    scores = np.zeros(4, dtype=float)
    counts = np.zeros(4, dtype=int)
    minima = np.zeros(4, dtype=int)
    for c in range(4):
        scores[c], counts[c], minima[c], _ = _score_one_class(
            pools[c], hams[c], p, top_m, radius
        )
    return _stable_argmax(scores), scores, counts, minima


def adaptive_soe_prediction(pools, hams, p, stages, confidence_threshold, min_ratio_threshold):
    """Adaptive logical-class mass estimation.

    Start with a cheap local stabilizer orbit. Expand only when the current top
    class is not sufficiently separated from the runner-up. Both posterior-like
    normalized confidence and top/runner-up ratio are reported; stopping requires
    either criterion to be strong enough.
    """
    history = []
    final = None
    for stage_idx, (top_m, radius) in enumerate(stages):
        pred, scores, counts, minima = fixed_soe_prediction(
            pools, hams, p, top_m=int(top_m), radius=int(radius)
        )
        order = np.argsort(scores)[::-1]
        top, second = float(scores[order[0]]), float(scores[order[1]])
        total = float(scores.sum())
        conf = top / total if total > 0 else 0.0
        ratio = top / max(second, 1e-300)
        final = (pred, scores, counts, minima, stage_idx, conf, ratio)
        history.append(
            {
                "stage": stage_idx,
                "top_m": int(top_m),
                "radius": int(radius),
                "prediction": pred,
                "confidence": conf,
                "top2_ratio": ratio,
                "expanded_total": int(counts.sum()),
            }
        )
        if stage_idx == len(stages) - 1:
            break
        if conf >= float(confidence_threshold) or ratio >= float(min_ratio_threshold):
            break
    return (*final, history)


def _cp_map_prediction(pools, hams):
    """Select the best minimum-weight configuration found across logical classes."""
    best = None
    for c in range(4):
        phys = hams[c].physical_from_bits(pools[c])
        n = hams[c].n_data
        wt = np.count_nonzero(phys[:, :n] | phys[:, n:], axis=1)
        w = int(wt.min())
        if best is None or w < best[0]:
            best = (w, c)
    return int(best[1]), int(best[0])


def _build_syndrome_models(code, oracle, cfg):
    """Precompute CP-SB pools/predictions for all d=3 syndrome keys."""
    models = {}
    rows = []
    for key in sorted(oracle.class_weight_counts):
        sx, sz = syndrome_from_key(key, code.Hz.shape[0], code.Hx.shape[0])
        pools, hams, ex0, ez0, stats = collect_cp_class_pools(
            code, oracle, sx, sz, cfg.solver_seed + 1000000 + key * 1009, cfg
        )
        cp_map_cls, cp_map_weight = _cp_map_prediction(pools, hams)
        fixed1_cls, fixed1_scores, fixed1_counts, _ = fixed_soe_prediction(pools, hams, cfg.p, 1, 1)
        fixed2_cls, fixed2_scores, fixed2_counts, _ = fixed_soe_prediction(pools, hams, cfg.p, 2, 2)
        (
            adaptive_cls,
            adaptive_scores,
            adaptive_counts,
            _,
            stage_idx,
            conf,
            top2_ratio,
            history,
        ) = adaptive_soe_prediction(
            pools,
            hams,
            cfg.p,
            cfg.stages,
            cfg.confidence_threshold,
            cfg.min_ratio_threshold,
        )
        exact_log, exact_masses = oracle.exact_mld_class(key, cfg.p)
        exact_cls = (int(exact_log[0]) << 1) | int(exact_log[1])
        exact_min = [int(np.flatnonzero(oracle.class_weight_counts[key][c])[0]) for c in range(4)]
        found_min = [stats[c]["best_weight"] for c in range(4)]
        models[key] = {
            "pools": pools,
            "hams": hams,
            "ex0": ex0,
            "ez0": ez0,
            "cp_map_cls": cp_map_cls,
            "fixed_r1_cls": fixed1_cls,
            "fixed_r2_cls": fixed2_cls,
            "adaptive_cls": adaptive_cls,
            "exact_cls": exact_cls,
            "adaptive_stage": int(stage_idx),
            "adaptive_conf": float(conf),
            "adaptive_top2_ratio": float(top2_ratio),
            "adaptive_counts": adaptive_counts,
            "history": history,
        }
        rows.append(
            {
                "syndrome_key": key,
                "cp_unique_candidates_total": sum(stats[c]["unique_cp_candidates"] for c in range(4)),
                "class_minima_all_exact": int(found_min == exact_min),
                "cp_map_class": cp_map_cls,
                "exact_mld_class": exact_cls,
                "fixed_r1_class": fixed1_cls,
                "fixed_r2_class": fixed2_cls,
                "adaptive_class": adaptive_cls,
                "adaptive_agrees_exact": int(adaptive_cls == exact_cls),
                "fixed_r1_agrees_exact": int(fixed1_cls == exact_cls),
                "fixed_r2_agrees_exact": int(fixed2_cls == exact_cls),
                "adaptive_stage": int(stage_idx),
                "adaptive_confidence": float(conf),
                "adaptive_top2_ratio": float(top2_ratio),
                "adaptive_expanded_total": int(adaptive_counts.sum()),
                "fixed_r1_expanded_total": int(fixed1_counts.sum()),
                "fixed_r2_expanded_total": int(fixed2_counts.sum()),
                "exact_winning_mass": float(exact_masses[exact_cls]),
                "adaptive_winning_score": float(adaptive_scores[adaptive_cls]),
            }
        )
    return models, pd.DataFrame(rows)


def run_cp_adaptive_confirmation(outdir, cfg=CPAdaptiveConfig()):
    from .oracle import build_d3_oracle, independent_decode, syndrome_key
    from .benchmark import wilson

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if cfg.distance != 3:
        raise ValueError("Exact-oracle CP-adaptive confirmation is intentionally d=3.")

    t0 = time.time()
    code = build_rotated_surface_code(cfg.distance)
    oracle = build_d3_oracle(code)
    models, syndrome_df = _build_syndrome_models(code, oracle, cfg)
    samples = sample_depolarizing(code, cfg.p, cfg.shots, cfg.seed)

    methods = [
        "matching",
        "single_MAP_config",
        "CP_SB_MAP",
        "CP_fixed_SOE_r1",
        "CP_fixed_SOE_r2",
        "CP_adaptive_SOE",
        "exact_logical_MLD",
    ]
    correct = {m: np.zeros(cfg.shots, dtype=np.uint8) for m in methods}
    shot_rows = []
    expanded = np.zeros(cfg.shots, dtype=int)
    stages = np.zeros(cfg.shots, dtype=int)

    for i in range(cfg.shots):
        sx, sz = samples["sx"][i], samples["sz"][i]
        key = syndrome_key(sx, sz)
        truth = int(samples["logical"][i, 0]) * 2 + int(samples["logical"][i, 1])
        model = models[key]

        ex0, ez0, match_log = independent_decode(code, oracle, sx, sz)
        match_cls = int(match_log[0]) * 2 + int(match_log[1])
        map_ex, map_ez = oracle.map_config[key]
        map_cls = _logical_class(code, map_ex, map_ez)
        pred = {
            "matching": match_cls,
            "single_MAP_config": map_cls,
            "CP_SB_MAP": model["cp_map_cls"],
            "CP_fixed_SOE_r1": model["fixed_r1_cls"],
            "CP_fixed_SOE_r2": model["fixed_r2_cls"],
            "CP_adaptive_SOE": model["adaptive_cls"],
            "exact_logical_MLD": model["exact_cls"],
        }
        for name in methods:
            correct[name][i] = int(pred[name] == truth)
        expanded[i] = int(model["adaptive_counts"].sum())
        stages[i] = int(model["adaptive_stage"])
        shot_rows.append(
            {
                "shot": i,
                "syndrome_key": key,
                "truth_class": truth,
                **{f"pred_{k}": v for k, v in pred.items()},
                **{f"correct_{k}": int(pred[k] == truth) for k in methods},
                "adaptive_stage": int(model["adaptive_stage"]),
                "adaptive_expanded_total": int(model["adaptive_counts"].sum()),
                "adaptive_confidence": float(model["adaptive_conf"]),
                "adaptive_top2_ratio": float(model["adaptive_top2_ratio"]),
            }
        )

    summary_rows = []
    for name in methods:
        k = int(correct[name].sum())
        lo, hi = wilson(k, cfg.shots)
        summary_rows.append(
            {
                "method": name,
                "logical_success": k / cfg.shots,
                "logical_error": 1.0 - k / cfg.shots,
                "ci95_lo": lo,
                "ci95_hi": hi,
            }
        )
    summary = pd.DataFrame(summary_rows)

    comparisons = []
    for a, b in [
        ("CP_adaptive_SOE", "matching"),
        ("CP_adaptive_SOE", "single_MAP_config"),
        ("CP_adaptive_SOE", "CP_SB_MAP"),
        ("CP_adaptive_SOE", "CP_fixed_SOE_r1"),
        ("CP_adaptive_SOE", "CP_fixed_SOE_r2"),
        ("exact_logical_MLD", "CP_adaptive_SOE"),
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
    comparisons = pd.DataFrame(comparisons)

    pd.DataFrame(shot_rows).to_csv(outdir / "shot_level.csv", index=False)
    summary.to_csv(outdir / "summary.csv", index=False)
    comparisons.to_csv(outdir / "paired_comparisons.csv", index=False)
    syndrome_df.to_csv(outdir / "syndrome_level.csv", index=False)

    meta = asdict(cfg)
    meta.update(
        {
            "elapsed_s": time.time() - t0,
            "all_256_syndromes_class_minima_exact_rate": float(syndrome_df["class_minima_all_exact"].mean()),
            "adaptive_exact_class_agreement_unweighted_syndromes": float(syndrome_df["adaptive_agrees_exact"].mean()),
            "fixed_r1_exact_class_agreement_unweighted_syndromes": float(syndrome_df["fixed_r1_agrees_exact"].mean()),
            "fixed_r2_exact_class_agreement_unweighted_syndromes": float(syndrome_df["fixed_r2_agrees_exact"].mean()),
            "adaptive_exact_class_agreement_shot_weighted": float(np.mean([
                models[syndrome_key(samples["sx"][i], samples["sz"][i])]["adaptive_cls"] ==
                models[syndrome_key(samples["sx"][i], samples["sz"][i])]["exact_cls"]
                for i in range(cfg.shots)
            ])),
            "mean_adaptive_expanded_candidates_per_shot": float(expanded.mean()),
            "median_adaptive_expanded_candidates_per_shot": float(np.median(expanded)),
            "fraction_stopped_after_stage0": float(np.mean(stages == 0)),
            "fraction_reached_stage1": float(np.mean(stages >= 1)),
            "fraction_reached_stage2": float(np.mean(stages >= 2)),
        }
    )
    (outdir / "metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    return summary, comparisons, syndrome_df, meta
