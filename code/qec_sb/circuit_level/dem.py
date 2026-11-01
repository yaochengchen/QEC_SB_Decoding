"""Detector error model -> matrices.

H  : detectors  x mechanisms (GF(2))         -- syndrome map
L  : observables x mechanisms (GF(2))        -- logical map
w  : log((1-p_i)/p_i) per mechanism          -- energy of mechanism i being "on"
tick: representative circuit tick per mechanism (for locality windows)
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp


@dataclass
class DEMMatrices:
    H: sp.csr_matrix          # (n_det, n_mech) uint8
    L: sp.csr_matrix          # (n_obs, n_mech) uint8
    probs: np.ndarray         # (n_mech,)
    weights: np.ndarray       # (n_mech,) float64, log((1-p)/p)
    tick: np.ndarray          # (n_mech,) int, representative tick (-1 if unknown)
    dem: object                # stim.DetectorErrorModel

    @property
    def n_mech(self):
        return self.H.shape[1]

    @property
    def n_det(self):
        return self.H.shape[0]

    @property
    def n_obs(self):
        return self.L.shape[0]

    def HL(self) -> sp.csr_matrix:
        return sp.vstack([self.H, self.L]).tocsr()


def dem_matrices(circuit: stim.Circuit, with_ticks: bool = True) -> DEMMatrices:
    """Undecomposed DEM (hyperedges kept) -> sparse matrices. Identical (D,L) sets are already
    merged by Stim; we additionally merge any residual duplicates by combining probabilities."""
    dem = circuit.detector_error_model(decompose_errors=False, flatten_loops=True,
                                       allow_gauge_detectors=False)
    n_det, n_obs = dem.num_detectors, dem.num_observables
    rows_d, rows_l, probs = [], [], []
    seen: dict[tuple, int] = {}
    for inst in dem.flattened():
        if inst.type != "error":
            continue
        p = inst.args_copy()[0]
        dets, obs = [], []
        for t in inst.targets_copy():
            if t.is_relative_detector_id():
                dets.append(t.val)
            elif t.is_logical_observable_id():
                obs.append(t.val)
        key = (tuple(sorted(set(_odd(dets)))), tuple(sorted(set(_odd(obs)))))
        if key in seen:                      # merge: p = p1(1-p2) + p2(1-p1)
            k = seen[key]
            probs[k] = probs[k] * (1 - p) + p * (1 - probs[k])
            continue
        seen[key] = len(probs)
        rows_d.append(key[0])
        rows_l.append(key[1])
        probs.append(p)
    n_mech = len(probs)
    H = _csr(rows_d, n_det, n_mech)
    L = _csr(rows_l, n_obs, n_mech)
    probs = np.asarray(probs)
    weights = np.log((1 - probs) / probs)
    tick = -np.ones(n_mech, dtype=np.int64)
    if with_ticks:
        tick = _ticks(circuit, dem, seen, n_mech)
    return DEMMatrices(H, L, probs, weights, tick, dem)


def _odd(ids):
    """Keep ids that appear an odd number of times (Stim can repeat a target)."""
    from collections import Counter
    return [i for i, c in Counter(ids).items() if c % 2]


def _csr(rows, nrows, ncols):
    indptr = [0]
    indices = []
    for cols in rows:
        indices.extend(cols)
        indptr.append(len(indices))
    # rows are mechanisms here; build (n_mech, nrows) then transpose
    M = sp.csr_matrix((np.ones(len(indices), dtype=np.uint8), indices, indptr), shape=(ncols, nrows))
    return M.T.tocsr()


def _ticks(circuit, dem, seen, n_mech):
    """Representative tick per mechanism using Stim's error explanations."""
    tick = -np.ones(n_mech, dtype=np.int64)
    try:
        expl = circuit.explain_detector_error_model_errors(dem_filter=dem, reduce_to_one_representative_error=True)
    except Exception:                      # old Stim without this API
        return tick
    for e in expl:
        dets, obs = [], []
        for t in e.dem_error_terms:
            d = t.dem_target
            if d.is_relative_detector_id():
                dets.append(d.val)
            elif d.is_logical_observable_id():
                obs.append(d.val)
        key = (tuple(sorted(set(_odd(dets)))), tuple(sorted(set(_odd(obs)))))
        k = seen.get(key)
        if k is None or not e.circuit_error_locations:
            continue
        tick[k] = min(loc.tick_offset for loc in e.circuit_error_locations)
    return tick


def sample(circuit: stim.Circuit, shots: int, seed: int):
    """Detector events (shots, n_det) bool and observable flips (shots, n_obs) bool."""
    s = circuit.compile_detector_sampler(seed=seed)
    det, obs = s.sample(shots=shots, separate_observables=True)
    return det, obs
