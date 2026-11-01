from __future__ import annotations

import math
import hashlib

import numpy as np
import pymatching

from ..rotated_code import build_rotated_surface_code, sample_depolarizing
from ..scalable import CSSMWPMDecoder, correction_class
from .statistics import exact_paired_statistics, wilson_interval as wilson


class PyMatchingCSSDecoder:
    """Official PyMatching decoder for the two CSS components.

    Under depolarizing code-capacity noise, each binary X or Z component has
    marginal error probability 2p/3. Correlations between the components are
    intentionally not supplied, matching the usual independent-CSS baseline.
    """

    def __init__(self, code, p: float):
        q = 2.0 * float(p) / 3.0
        probabilities = np.full(code.n, q, dtype=float)
        self.code = code
        self.x_matching = pymatching.Matching(code.Hz, error_probabilities=probabilities)
        self.z_matching = pymatching.Matching(code.Hx, error_probabilities=probabilities)

    def decode(self, sx, sz):
        ex = np.asarray(self.x_matching.decode(np.asarray(sx, dtype=np.uint8)), dtype=np.uint8)
        ez = np.asarray(self.z_matching.decode(np.asarray(sz, dtype=np.uint8)), dtype=np.uint8)
        return ex, ez


def _bp_matching_component(H, syndrome, q, iterations, damping):
    H = np.asarray(H, dtype=np.uint8)
    syndrome = np.asarray(syndrome, dtype=np.uint8)
    m, n = H.shape
    prior = math.log((1.0 - q) / q)
    v_to_c = np.zeros((m, n), dtype=float)
    c_to_v = np.zeros((m, n), dtype=float)
    v_to_c[H == 1] = prior
    checks = [np.flatnonzero(H[i]).tolist() for i in range(m)]
    variables = [np.flatnonzero(H[:, j]).tolist() for j in range(n)]

    for _ in range(int(iterations)):
        next_c = np.zeros_like(c_to_v)
        for i, neighbors in enumerate(checks):
            sign = -1.0 if syndrome[i] else 1.0
            for j in neighbors:
                product = sign
                for other in neighbors:
                    if other != j:
                        product *= math.tanh(v_to_c[i, other] / 2.0)
                product = float(np.clip(product, -1.0 + 1e-12, 1.0 - 1e-12))
                next_c[i, j] = 2.0 * math.atanh(product)
        c_to_v = damping * c_to_v + (1.0 - damping) * next_c
        next_v = np.zeros_like(v_to_c)
        for j, neighbors in enumerate(variables):
            for i in neighbors:
                next_v[i, j] = prior + sum(c_to_v[other, j] for other in neighbors if other != i)
        v_to_c = damping * v_to_c + (1.0 - damping) * next_v

    posterior = np.full(n, prior, dtype=float)
    for j, neighbors in enumerate(variables):
        posterior[j] += sum(c_to_v[i, j] for i in neighbors)
    estimate = (posterior < 0).astype(np.uint8)
    residual = syndrome ^ ((H @ estimate) & 1)
    if residual.any():
        residual_matching = pymatching.Matching(H, weights=np.maximum(np.abs(posterior), 1e-6))
        estimate ^= np.asarray(residual_matching.decode(residual), dtype=np.uint8)
    if not np.array_equal((H @ estimate) & 1, syndrome):
        raise RuntimeError("BP+Matching produced an incompatible correction")
    return estimate


class BPMatchingCSSDecoder:
    """Fixed-parameter sum-product BP followed by residual weighted matching."""

    def __init__(self, code, p: float, iterations: int = 12, damping: float = 0.4):
        self.code = code
        self.q = 2.0 * float(p) / 3.0
        self.iterations = int(iterations)
        self.damping = float(damping)

    def decode(self, sx, sz):
        ex = _bp_matching_component(self.code.Hz, sx, self.q, self.iterations, self.damping)
        ez = _bp_matching_component(self.code.Hx, sz, self.q, self.iterations, self.damping)
        return ex, ez


def evaluate_code_capacity_baselines(distance, p, shots, seed, bp_iterations=12, bp_damping=0.4):
    code = build_rotated_surface_code(int(distance))
    samples = sample_depolarizing(code, float(p), int(shots), int(seed))
    decoders = {
        "NetworkX_MWPM": CSSMWPMDecoder(code),
        "PyMatching": PyMatchingCSSDecoder(code, p),
        "BP_Matching": BPMatchingCSSDecoder(code, p, bp_iterations, bp_damping),
    }
    correct = {name: [] for name in decoders}
    predictions = {name: [] for name in decoders}
    truth_classes = []
    syndrome_keys = []
    for index in range(int(shots)):
        truth = int(samples["logical"][index, 0]) * 2 + int(samples["logical"][index, 1])
        sx, sz = samples["sx"][index], samples["sz"][index]
        truth_classes.append(truth)
        key = 0
        for bit in np.concatenate([sx, sz]):
            key = (key << 1) | int(bit)
        syndrome_keys.append(key)
        for name, decoder in decoders.items():
            ex, ez = decoder.decode(sx, sz)
            prediction = int(correction_class(code, ex, ez))
            predictions[name].append(prediction)
            correct[name].append(int(prediction == truth))
    digest = hashlib.sha256()
    for key in ("pauli", "sx", "sz", "logical"):
        digest.update(np.ascontiguousarray(samples[key]).tobytes())
    return {
        "distance": int(distance),
        "p": float(p),
        "shots": int(shots),
        "sample_seed": int(seed),
        "bp_iterations": int(bp_iterations),
        "bp_damping": float(bp_damping),
        "correct": correct,
        "predictions": predictions,
        "sample_sha256": digest.hexdigest(),
        "truth_classes": truth_classes,
        "syndrome_keys": syndrome_keys,
    }


def summarize_paired_methods(correct, reference="QEC_SB"):
    reference_values = np.asarray(correct[reference], dtype=np.uint8)
    rows = []
    for method, values in correct.items():
        values = np.asarray(values, dtype=np.uint8)
        if len(values) != len(reference_values):
            raise ValueError("paired method arrays must have equal length")
        paired = exact_paired_statistics(reference_values, values)
        delta = reference_values.astype(float) - values.astype(float)
        success = float(values.mean())
        lo, hi = wilson(int(values.sum()), len(values))
        rows.append(
            {
                "method": method,
                "shots": int(len(values)),
                "logical_success": success,
                "logical_error": 1.0 - success,
                "ci95_lo": lo,
                "ci95_hi": hi,
                "reference": reference,
                "reference_minus_method": float(delta.mean()),
                "paired_ci95_lo": paired["paired_ci95_lo"],
                "paired_ci95_hi": paired["paired_ci95_hi"],
                "reference_correct_method_wrong": int(np.sum((reference_values == 1) & (values == 0))),
                "reference_wrong_method_correct": int(np.sum((reference_values == 0) & (values == 1))),
                "discordant_wins": paired["discordant_wins"],
                "discordant_losses": paired["discordant_losses"],
                "exact_p_value": paired["exact_p_value"],
            }
        )
    return rows

