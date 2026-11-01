from __future__ import annotations

import numpy as np
from scipy.stats import beta, binomtest


def adjust_p_values(p_values, method: str) -> list[float]:
    """Adjust a family of p-values using Holm or Bonferroni correction."""
    if method not in {"holm", "bonferroni"}:
        raise ValueError("method must be 'holm' or 'bonferroni'")
    try:
        if np.iscomplexobj(p_values):
            raise ValueError("p_values must contain real numbers")
        values = np.asarray(p_values, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("p_values must be a one-dimensional sequence of finite numbers in [0, 1]") from exc
    if values.ndim != 1 or not np.all(np.isfinite(values)) or np.any((values < 0.0) | (values > 1.0)):
        raise ValueError("p_values must be a one-dimensional sequence of finite numbers in [0, 1]")

    count = len(values)
    if method == "bonferroni":
        adjusted = np.minimum(count * values, 1.0)
    else:
        order = np.argsort(values, kind="stable")
        sorted_values = values[order]
        multipliers = count - np.arange(count)
        sorted_adjusted = np.minimum(np.maximum.accumulate(sorted_values * multipliers), 1.0)
        adjusted = np.empty(count, dtype=float)
        adjusted[order] = sorted_adjusted
    return [float(value) for value in adjusted]


def wilson_interval(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    denominator = 1.0 + z * z / n
    center = (p + z * z / (2.0 * n)) / denominator
    half = z * np.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n)) / denominator
    return float(center - half), float(center + half)


def _clopper_pearson(k: int, n: int, alpha: float) -> tuple[float, float]:
    if n <= 0:
        raise ValueError("n must be positive")
    lo = 0.0 if k == 0 else float(beta.ppf(alpha / 2.0, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1.0 - alpha / 2.0, k + 1, n - k))
    return lo, hi


def exact_paired_statistics(method, baseline, confidence: float = 0.95) -> dict:
    """Exact McNemar/sign inference for paired binary outcomes.

    The difference interval is a conservative simultaneous interval: separate
    Clopper-Pearson intervals for win and loss proportions, with Bonferroni
    allocation across the two proportions, are differenced.
    """
    method = np.asarray(method, dtype=np.uint8)
    baseline = np.asarray(baseline, dtype=np.uint8)
    if method.ndim != 1 or baseline.ndim != 1 or len(method) != len(baseline) or len(method) == 0:
        raise ValueError("paired samples must be one-dimensional and have the same positive length")
    if np.any(method > 1) or np.any(baseline > 1):
        raise ValueError("paired outcomes must be binary")
    n = len(method)
    wins = int(np.sum((method == 1) & (baseline == 0)))
    losses = int(np.sum((method == 0) & (baseline == 1)))
    discordant = wins + losses
    p_value = 1.0 if discordant == 0 else float(binomtest(wins, discordant, 0.5).pvalue)
    family_alpha = 1.0 - float(confidence)
    # Each interval has coverage 1-alpha/2, hence joint coverage >=1-alpha.
    win_lo, win_hi = _clopper_pearson(wins, n, family_alpha / 2.0)
    loss_lo, loss_hi = _clopper_pearson(losses, n, family_alpha / 2.0)
    return {
        "shots": int(n),
        "method_success": float(method.mean()),
        "baseline_success": float(baseline.mean()),
        "paired_difference": float((method.astype(float) - baseline.astype(float)).mean()),
        "paired_ci95_lo": float(win_lo - loss_hi),
        "paired_ci95_hi": float(win_hi - loss_lo),
        "discordant_wins": wins,
        "discordant_losses": losses,
        "exact_p_value": p_value,
    }
