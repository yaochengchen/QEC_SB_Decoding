from dataclasses import dataclass, asdict
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from .rotated_code import build_rotated_surface_code, sample_depolarizing
from .oracle import build_d3_oracle, independent_decode, syndrome_key
from .p3 import _logical_class, _paired_diff_ci, _stable_argmax
from .benchmark import wilson
from .adaptive_soe import (
    CPAdaptiveConfig,
    syndrome_from_key,
    collect_cp_class_pools,
    fixed_soe_prediction,
)
from .cp_sb import stabilizer_hamming_neighborhood


@dataclass
class ImportanceMassConfig:
    """P4.1 adaptive logical-mass estimator settings.

    Full stabilizer-orbit enumeration is forbidden in the estimator.  A cheap
    radius-1 SOE gate handles easy syndromes.  Ambiguous syndromes are passed to
    a stratified multiple-importance sampler (MIS): a small deterministic local
    core is summed exactly and the remaining tail mass is estimated unbiasedly.
    """

    # Cheap stage inherited from P3.3.
    stage0_top_m: int = 1
    stage0_radius: int = 1
    stage0_ratio_threshold: float = 1.05

    # Deterministic control-variate/core around multiple CP-SB low-energy seeds.
    top_centers_per_class: int = 2
    core_radius: int = 2

    # Shared mixture proposal over stabilizer coordinates u.
    proposal_rhos: tuple = (0.08, 0.20, 0.40, 0.50)
    uniform_mix: float = 0.05

    # Sequential paired-CI stopping.
    batch_size: int = 32
    min_samples: int = 64
    max_samples: int = 160
    z_score: float = 2.58

    # Hard P4.1 requirement: never deliberately cover a full logical orbit.
    prohibit_full_enumeration: bool = True


@dataclass
class P41Config:
    distance: int = 3
    p: float = 0.08
    shots: int = 10000
    seed: int = 20261018
    solver_seed: int = 20261014
    replicas: int = 8
    steps: int = 120
    dt: float = 0.45
    coupling_scale: float = 0.70
    warm_strength: float = 0.08
    importance: ImportanceMassConfig = None

    def __post_init__(self):
        if self.importance is None:
            self.importance = ImportanceMassConfig()


def _encode_bits(u):
    """Return hashable state keys without truncating any stabilizer bits."""
    u = np.asarray(u, dtype=np.uint8)
    one = u.ndim == 1
    uu = u[None, :] if one else u
    ids = [row.tobytes() for row in np.packbits(uu, axis=1, bitorder="little")]
    return ids[0] if one else ids


def _pauli_weights(ham, u):
    physical = ham.physical_from_bits(u)
    one = physical.ndim == 1
    pp = physical[None, :] if one else physical
    n = ham.n_data
    w = np.count_nonzero(pp[:, :n] | pp[:, n:], axis=1)
    return int(w[0]) if one else w


def _best_centers(pool_u, ham, k):
    w = _pauli_weights(ham, pool_u)
    order = np.argsort(w, kind="stable")[: max(1, min(int(k), len(pool_u)))]
    return pool_u[order].copy()


def _build_core_and_proposal_centers(pools, hams, p, cfg: ImportanceMassConfig):
    ratio = (p / 3.0) / (1.0 - p)
    class_cores = []
    class_core_ids = []
    core_masses = np.zeros(4, dtype=float)
    all_centers = []

    for c in range(4):
        centers = _best_centers(pools[c], hams[c], cfg.top_centers_per_class)
        all_centers.extend(centers)
        neighborhoods = [
            stabilizer_hamming_neighborhood(center, hams[c].n_spins, cfg.core_radius)
            for center in centers
        ]
        core = np.unique(np.vstack(neighborhoods).astype(np.uint8), axis=0)
        class_cores.append(core)
        class_core_ids.append(set(_encode_bits(core)))
        w = _pauli_weights(hams[c], core)
        core_masses[c] = float(np.sum(ratio ** w))

    proposal_centers = np.unique(np.asarray(all_centers, dtype=np.uint8), axis=0)
    return class_cores, class_core_ids, core_masses, proposal_centers


def _proposal_probability(u, centers, rhos, uniform_mix):
    u = np.asarray(u, dtype=np.uint8)
    one = u.ndim == 1
    uu = u[None, :] if one else u
    m = uu.shape[1]
    alpha = float(uniform_mix)
    if not (0.0 < alpha < 1.0):
        raise ValueError("uniform_mix must lie strictly between 0 and 1.")
    rhos = np.asarray(rhos, dtype=float)
    if np.any((rhos <= 0.0) | (rhos >= 0.5 + 1e-15)):
        raise ValueError("proposal_rhos must be in (0, 0.5].")

    g = np.full(len(uu), alpha / float(2 ** m), dtype=float)
    local = np.zeros(len(uu), dtype=float)
    for center in centers:
        h = np.count_nonzero(uu ^ center[None, :], axis=1)
        for rho in rhos:
            local += (rho ** h) * ((1.0 - rho) ** (m - h))
    local /= float(len(centers) * len(rhos))
    g += (1.0 - alpha) * local
    return float(g[0]) if one else g


def _sample_proposal(rng, n_samples, n_spins, centers, rhos, uniform_mix):
    n_samples = int(n_samples)
    u = np.empty((n_samples, n_spins), dtype=np.uint8)
    alpha = float(uniform_mix)
    choose_uniform = rng.random(n_samples) < alpha
    n_uniform = int(choose_uniform.sum())
    if n_uniform:
        u[choose_uniform] = rng.integers(
            0, 2, size=(n_uniform, n_spins), dtype=np.uint8
        )
    ids = np.flatnonzero(~choose_uniform)
    if len(ids):
        cids = rng.integers(0, len(centers), size=len(ids))
        rids = rng.integers(0, len(rhos), size=len(ids))
        base = centers[cids].copy()
        probs = np.asarray(rhos, dtype=float)[rids][:, None]
        flips = (rng.random((len(ids), n_spins)) < probs).astype(np.uint8)
        u[ids] = base ^ flips
    return u


def _paired_lower_bounds(core_masses, tail_contrib, z_score):
    """Return winner and paired lower bounds for Z_winner - Z_competitor."""
    estimates = core_masses + tail_contrib.mean(axis=0)
    top = int(np.argmax(estimates))
    lowers = np.full(4, np.inf, dtype=float)
    deltas = np.zeros(4, dtype=float)
    ses = np.zeros(4, dtype=float)
    n = len(tail_contrib)
    for c in range(4):
        if c == top:
            continue
        diff_samples = tail_contrib[:, top] - tail_contrib[:, c]
        delta = float(core_masses[top] - core_masses[c] + diff_samples.mean())
        se = float(diff_samples.std(ddof=1) / np.sqrt(n)) if n > 1 else np.inf
        lowers[c] = delta - float(z_score) * se
        deltas[c] = delta
        ses[c] = se
    return top, estimates, lowers, deltas, ses


def importance_mass_prediction(pools, hams, p, rng_seed, cfg=ImportanceMassConfig()):
    """Estimate four logical-class masses without full orbit enumeration.

    The estimator is unbiased for every Z_L after adding the exact deterministic
    local core: samples u ~ g(u) are shared by all logical classes and contribute

        1[u not in core_L] * r**weight_L(u) / g(u).

    Sharing u gives a paired estimator of Z_a-Z_b, which is what decoding needs.
    The final ambiguity rule chooses the smallest class among statistically active
    classes; exact ties are physically equivalent MLD decisions.
    """
    n_spins = hams[0].n_spins
    orbit_size = 2 ** n_spins
    class_cores, class_core_ids, core_masses, centers = _build_core_and_proposal_centers(
        pools, hams, p, cfg
    )
    max_core = max(len(v) for v in class_cores)
    if cfg.prohibit_full_enumeration and cfg.max_samples + max_core >= orbit_size:
        raise ValueError(
            "P4.1 full-enumeration guard failed: max_samples + largest core "
            f"must be < orbit size ({cfg.max_samples}+{max_core} >= {orbit_size})."
        )

    rng = np.random.default_rng(int(rng_seed))
    ratio = (p / 3.0) / (1.0 - p)
    batches = []
    unique_sample_ids = set()
    stop_reason = "max_samples_ambiguous"
    final = None

    for start in range(0, int(cfg.max_samples), int(cfg.batch_size)):
        n_batch = min(int(cfg.batch_size), int(cfg.max_samples) - start)
        u = _sample_proposal(
            rng,
            n_batch,
            n_spins,
            centers,
            cfg.proposal_rhos,
            cfg.uniform_mix,
        )
        sample_ids = _encode_bits(u)
        unique_sample_ids.update(sample_ids)
        g = _proposal_probability(u, centers, cfg.proposal_rhos, cfg.uniform_mix)
        y = np.zeros((len(u), 4), dtype=float)
        for c in range(4):
            w = _pauli_weights(hams[c], u)
            contrib = (ratio ** w) / g
            in_core = np.fromiter(
                (v in class_core_ids[c] for v in sample_ids),
                dtype=bool,
                count=len(sample_ids),
            )
            contrib[in_core] = 0.0
            y[:, c] = contrib
        batches.append(y)
        tail = np.vstack(batches)

        if len(tail) < int(cfg.min_samples):
            continue
        top, estimates, lowers, deltas, ses = _paired_lower_bounds(
            core_masses, tail, cfg.z_score
        )
        resolved = all(lowers[c] > 0.0 for c in range(4) if c != top)
        final = (top, estimates, lowers, deltas, ses, tail)
        if resolved:
            stop_reason = "paired_ci_resolved"
            break

    if final is None:
        tail = np.vstack(batches)
        final = (*_paired_lower_bounds(core_masses, tail, cfg.z_score), tail)

    top, estimates, lowers, deltas, ses, tail = final
    active = [top]
    if stop_reason != "paired_ci_resolved":
        for c in range(4):
            if c == top:
                continue
            if lowers[c] <= 0.0:
                active.append(c)
        prediction = int(min(active))
    else:
        prediction = int(top)

    max_seen_per_class = max_core + len(unique_sample_ids)
    if cfg.prohibit_full_enumeration and max_seen_per_class >= orbit_size:
        raise RuntimeError("Estimator accidentally covered a full logical orbit.")

    diagnostics = {
        "prediction": prediction,
        "raw_top_class": int(top),
        "active_classes": tuple(int(v) for v in sorted(set(active))),
        "mass_estimates": estimates,
        "core_masses": core_masses,
        "paired_lower_bounds": lowers,
        "paired_deltas": deltas,
        "paired_ses": ses,
        "samples_drawn": int(len(tail)),
        "unique_shared_samples": int(len(unique_sample_ids)),
        "core_candidates_total": int(sum(len(v) for v in class_cores)),
        "max_core_candidates_one_class": int(max_core),
        "proposal_centers": int(len(centers)),
        "stop_reason": stop_reason,
        "orbit_size_per_class": int(orbit_size),
        "full_orbit_evaluations": int(4 * orbit_size),
        "candidate_evaluations": int(sum(len(v) for v in class_cores) + 4 * len(tail)),
        "max_seen_states_one_class_upper_bound": int(max_seen_per_class),
    }
    return prediction, diagnostics


def adaptive_p41_prediction(pools, hams, p, rng_seed, cfg=ImportanceMassConfig()):
    """Cheap radius-1 gate, then importance sampling only for ambiguous syndromes."""
    stage0_cls, scores, counts, minima = fixed_soe_prediction(
        pools,
        hams,
        p,
        top_m=cfg.stage0_top_m,
        radius=cfg.stage0_radius,
    )
    order = np.argsort(scores)[::-1]
    top = float(scores[order[0]])
    second = float(scores[order[1]])
    ratio = top / max(second, 1e-300)
    stage0_cost = int(counts.sum())
    if ratio >= float(cfg.stage0_ratio_threshold):
        return int(stage0_cls), {
            "stage": 0,
            "stage0_prediction": int(stage0_cls),
            "stage0_top2_ratio": float(ratio),
            "samples_drawn": 0,
            "unique_shared_samples": 0,
            "core_candidates_total": stage0_cost,
            "candidate_evaluations": stage0_cost,
            "full_orbit_evaluations": int(4 * (2 ** hams[0].n_spins)),
            "stop_reason": "stage0_ratio_resolved",
            "active_classes": (int(stage0_cls),),
            "orbit_size_per_class": int(2 ** hams[0].n_spins),
            "max_seen_states_one_class_upper_bound": int(max(counts)),
        }

    pred, diag = importance_mass_prediction(pools, hams, p, rng_seed, cfg)
    diag.update(
        {
            "stage": 1,
            "stage0_prediction": int(stage0_cls),
            "stage0_top2_ratio": float(ratio),
        }
    )
    return pred, diag


def _exact_optimal_set(oracle, key, p):
    counts = oracle.class_weight_counts[key]
    ratio = (p / 3.0) / (1.0 - p)
    masses = counts @ (ratio ** np.arange(counts.shape[1]))
    mx = float(masses.max())
    optimal = tuple(
        int(v) for v in np.flatnonzero(np.isclose(masses, mx, rtol=1e-12, atol=1e-15))
    )
    return optimal, masses


def _build_p41_models(code, oracle, cfg: P41Config):
    cp_cfg = CPAdaptiveConfig(
        distance=cfg.distance,
        p=cfg.p,
        shots=cfg.shots,
        seed=cfg.seed,
        solver_seed=cfg.solver_seed,
        replicas=cfg.replicas,
        steps=cfg.steps,
        dt=cfg.dt,
        coupling_scale=cfg.coupling_scale,
        warm_strength=cfg.warm_strength,
    )
    models = {}
    rows = []
    for key in sorted(oracle.class_weight_counts):
        sx, sz = syndrome_from_key(key, code.Hz.shape[0], code.Hx.shape[0])
        pools, hams, ex0, ez0, stats = collect_cp_class_pools(
            code,
            oracle,
            sx,
            sz,
            cfg.solver_seed + 1000000 + key * 1009,
            cp_cfg,
        )
        pred, diag = adaptive_p41_prediction(
            pools,
            hams,
            cfg.p,
            rng_seed=cfg.seed + 7000000 + key * 7919,
            cfg=cfg.importance,
        )
        exact_log, exact_masses = oracle.exact_mld_class(key, cfg.p)
        exact_cls = (int(exact_log[0]) << 1) | int(exact_log[1])
        optimal_set, masses = _exact_optimal_set(oracle, key, cfg.p)
        unique_winner = len(optimal_set) == 1
        models[key] = {
            "prediction": int(pred),
            "diag": diag,
            "exact_cls": int(exact_cls),
            "optimal_set": optimal_set,
            "exact_masses": masses,
        }
        rows.append(
            {
                "syndrome_key": int(key),
                "prediction": int(pred),
                "exact_stable_class": int(exact_cls),
                "exact_optimal_classes": ";".join(str(v) for v in optimal_set),
                "n_exact_optimal_classes": len(optimal_set),
                "unique_exact_winner": int(unique_winner),
                "strict_exact_class_agreement": int(pred == exact_cls),
                "optimal_set_agreement": int(pred in optimal_set),
                "stage": int(diag["stage"]),
                "stage0_top2_ratio": float(diag["stage0_top2_ratio"]),
                "samples_drawn": int(diag["samples_drawn"]),
                "unique_shared_samples": int(diag["unique_shared_samples"]),
                "candidate_evaluations": int(diag["candidate_evaluations"]),
                "full_orbit_evaluations": int(diag["full_orbit_evaluations"]),
                "evaluation_fraction_of_full": float(
                    diag["candidate_evaluations"] / diag["full_orbit_evaluations"]
                ),
                "stop_reason": diag["stop_reason"],
                "active_classes": ";".join(str(v) for v in diag["active_classes"]),
            }
        )
    return models, pd.DataFrame(rows)


def run_p41_confirmation(outdir, cfg=P41Config()):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if cfg.distance != 3:
        raise ValueError("P4.1 oracle confirmation is intentionally d=3.")

    t0 = time.time()
    code = build_rotated_surface_code(cfg.distance)
    oracle = build_d3_oracle(code)
    models, syndrome_df = _build_p41_models(code, oracle, cfg)
    samples = sample_depolarizing(code, cfg.p, cfg.shots, cfg.seed)

    methods = ["matching", "single_MAP_config", "P4_1_importance", "exact_logical_MLD"]
    correct = {m: np.zeros(cfg.shots, dtype=np.uint8) for m in methods}
    shot_rows = []

    for i in range(cfg.shots):
        sx, sz = samples["sx"][i], samples["sz"][i]
        key = syndrome_key(sx, sz)
        truth = int(samples["logical"][i, 0]) * 2 + int(samples["logical"][i, 1])
        model = models[key]
        _, _, match_log = independent_decode(code, oracle, sx, sz)
        match_cls = int(match_log[0]) * 2 + int(match_log[1])
        map_ex, map_ez = oracle.map_config[key]
        map_cls = _logical_class(code, map_ex, map_ez)
        pred = {
            "matching": match_cls,
            "single_MAP_config": map_cls,
            "P4_1_importance": int(model["prediction"]),
            "exact_logical_MLD": int(model["exact_cls"]),
        }
        for name in methods:
            correct[name][i] = int(pred[name] == truth)
        diag = model["diag"]
        shot_rows.append(
            {
                "shot": i,
                "syndrome_key": int(key),
                "truth_class": int(truth),
                **{f"pred_{k}": int(v) for k, v in pred.items()},
                **{f"correct_{k}": int(pred[k] == truth) for k in methods},
                "prediction_is_exact_optimal": int(pred["P4_1_importance"] in model["optimal_set"]),
                "stage": int(diag["stage"]),
                "samples_drawn": int(diag["samples_drawn"]),
                "candidate_evaluations": int(diag["candidate_evaluations"]),
                "full_orbit_evaluations": int(diag["full_orbit_evaluations"]),
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
        ("P4_1_importance", "matching"),
        ("P4_1_importance", "single_MAP_config"),
        ("exact_logical_MLD", "P4_1_importance"),
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
    shot_df = pd.DataFrame(shot_rows)

    # Exact expected logical success under the known code-capacity depolarizing
    # model. This removes finite-shot fluctuations and arbitrary tie-breaking
    # effects from the decoder comparison. The oracle masses omit a common
    # factor (1-p)^n, restored here.
    common = (1.0 - cfg.p) ** code.n
    bayes_acc = {m: 0.0 for m in methods}
    probability_check = 0.0
    for key, model in models.items():
        sx, sz = syndrome_from_key(key, code.Hz.shape[0], code.Hx.shape[0])
        _, _, match_log = independent_decode(code, oracle, sx, sz)
        match_cls = int(match_log[0]) * 2 + int(match_log[1])
        map_ex, map_ez = oracle.map_config[key]
        map_cls = _logical_class(code, map_ex, map_ez)
        preds = {
            "matching": match_cls,
            "single_MAP_config": map_cls,
            "P4_1_importance": int(model["prediction"]),
            "exact_logical_MLD": int(model["exact_cls"]),
        }
        masses = np.asarray(model["exact_masses"], dtype=float)
        probability_check += common * float(masses.sum())
        for name, pred_cls in preds.items():
            bayes_acc[name] += common * float(masses[pred_cls])
    bayes_expected = pd.DataFrame([
        {
            "method": name,
            "expected_logical_success": float(bayes_acc[name]),
            "expected_logical_error": float(1.0 - bayes_acc[name]),
        }
        for name in methods
    ])

    summary.to_csv(outdir / "summary.csv", index=False)
    comparisons.to_csv(outdir / "paired_comparisons.csv", index=False)
    bayes_expected.to_csv(outdir / "bayes_expected_summary.csv", index=False)
    syndrome_df.to_csv(outdir / "syndrome_level.csv", index=False)
    shot_df.to_csv(outdir / "shot_level.csv", index=False)

    unique_df = syndrome_df[syndrome_df["unique_exact_winner"] == 1]
    meta = asdict(cfg)
    # dataclasses.asdict has already converted nested ImportanceMassConfig.
    meta.update(
        {
            "elapsed_s": time.time() - t0,
            "n_syndromes": int(len(syndrome_df)),
            "n_unique_winner_syndromes": int(len(unique_df)),
            "n_tied_optimum_syndromes": int(len(syndrome_df) - len(unique_df)),
            "strict_exact_class_agreement_unweighted": float(
                syndrome_df["strict_exact_class_agreement"].mean()
            ),
            "optimal_set_agreement_unweighted": float(
                syndrome_df["optimal_set_agreement"].mean()
            ),
            "unique_winner_strict_agreement": float(
                unique_df["strict_exact_class_agreement"].mean()
            ),
            "fraction_syndromes_using_importance": float(np.mean(syndrome_df["stage"] == 1)),
            "fraction_shots_using_importance": float(np.mean(shot_df["stage"] == 1)),
            "mean_samples_per_importance_syndrome": float(
                syndrome_df.loc[syndrome_df["stage"] == 1, "samples_drawn"].mean()
            ),
            "mean_candidate_evaluations_per_shot": float(shot_df["candidate_evaluations"].mean()),
            "median_candidate_evaluations_per_shot": float(shot_df["candidate_evaluations"].median()),
            "mean_fraction_of_full_orbit_cost_per_shot": float(
                (shot_df["candidate_evaluations"] / shot_df["full_orbit_evaluations"]).mean()
            ),
            "max_candidate_evaluations_any_syndrome": int(syndrome_df["candidate_evaluations"].max()),
            "full_orbit_evaluations_per_syndrome": int(syndrome_df["full_orbit_evaluations"].iloc[0]),
            "algorithm_used_full_enumeration": False,
            "probability_normalization_check": float(probability_check),
            "P4_1_expected_success": float(bayes_acc["P4_1_importance"]),
            "exact_MLD_expected_success": float(bayes_acc["exact_logical_MLD"]),
            "expected_success_regret_vs_exact": float(
                bayes_acc["exact_logical_MLD"] - bayes_acc["P4_1_importance"]
            ),
        }
    )
    (outdir / "metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    return summary, comparisons, syndrome_df, shot_df, meta
