"""Exact local stabilizer-orbit expansion (SOE) for radius 1 and 2 in O(m) per class.

Same quantity as ``adaptive_soe.fixed_soe_prediction(pools, hams, p, top_m=1, radius=r)``:

    Z_hat_L = sum_{u : d_H(u, u*) <= r} w^{|e_L + u G|},   w = (p/3)/(1-p),

where u* is the lowest-weight candidate of class L (first one under a stable sort, exactly as in
``adaptive_soe._score_one_class``). The old code enumerates the 1 + m + m(m-1)/2 states of the
ball and recomputes every Pauli weight (O(m^2 n)). Here:

  * The Pauli weight is a sum over qubits, so flipping generator j changes it by
    dE_j = sum_{q in supp j} ([new e_q != I] - [old e_q != I])  (only supp j is touched).
  * If supp j and supp k share no qubit, dE_{jk} = dE_j + dE_k exactly.
  * Hence the number of ball states at each weight is an integer histogram:
        center                      -> 1 state at dE = 0
        single flips                -> histogram of dE_j
        all pairs, as if disjoint   -> (h*h - diag)/2   (h = histogram of dE_j, * = convolution)
        overlapping pairs (O(m))    -> remove dE_j + dE_k, add the true dE_{jk}
    and Z_hat = sum_t N_t w^(E0 + t).

The histogram is exact integer arithmetic, so the result is the same number as the enumeration
(up to floating-point summation order, ~1e-15 relative). Classes that are tied exactly (identical
histograms) give bit-identical floats and are resolved to the lowest class index by
``_stable_argmax``, the convention the old code intends. Only the generator matrix G is used;
the pair list is built from G's qubit supports, so nothing here is specific to the surface code.
The per-code structure (supports, overlapping pairs) is built once per G and cached.
"""

from __future__ import annotations

import numpy as np

from .p3 import _stable_argmax

_CACHE: dict = {}


class SOEStructure:
    """Per-code precomputation: padded generator supports and the overlapping generator pairs."""

    def __init__(self, G: np.ndarray, n: int):
        G = np.asarray(G, dtype=np.uint8)
        m = int(G.shape[0])
        gx, gz = G[:, :n], G[:, n:]
        supp = (gx | gz).astype(bool)                       # [m, n] qubit support of each generator
        K = max(1, int(supp.sum(axis=1).max()))
        self.m, self.n, self.K = m, int(n), K
        self.G_f = G.astype(np.float64)                     # BLAS matmul; sums <= m are exact in float64

        # Single generators, padded to K qubits. Padding uses qubit 0 with zero flips -> contributes 0.
        self.q1 = np.zeros((m, K), dtype=np.int64)
        self.fx1 = np.zeros((m, K), dtype=np.uint8)
        self.fz1 = np.zeros((m, K), dtype=np.uint8)
        for j in range(m):
            qs = np.flatnonzero(supp[j])
            self.q1[j, :len(qs)] = qs
            self.fx1[j, :len(qs)] = gx[j, qs]
            self.fz1[j, :len(qs)] = gz[j, qs]

        # Pairs (j < k) whose qubit supports intersect: the only pairs where dE is not additive.
        ov = (supp.astype(np.int32) @ supp.T.astype(np.int32)) > 0
        pj, pk = np.nonzero(np.triu(ov, 1))
        self.pj, self.pk = pj.astype(np.int64), pk.astype(np.int64)
        U = 2 * K
        P = len(pj)
        self.q2 = np.zeros((P, U), dtype=np.int64)
        self.fx2 = np.zeros((P, U), dtype=np.uint8)
        self.fz2 = np.zeros((P, U), dtype=np.uint8)
        for t, (j, k) in enumerate(zip(pj, pk)):
            qs = np.flatnonzero(supp[j] | supp[k])
            self.q2[t, :len(qs)] = qs
            self.fx2[t, :len(qs)] = gx[j, qs] ^ gx[k, qs]
            self.fz2[t, :len(qs)] = gz[j, qs] ^ gz[k, qs]
        self.n_overlap_pairs = int(P)
        # histogram index t <-> weight change (t - OFF); |dE_jk| <= |supp j U supp k| <= 2K
        self.OFF = 2 * K
        self.L = 4 * K + 1


def get_structure(ham) -> SOEStructure:
    G = ham.generators
    key = (int(ham.n_data), G.shape, G.tobytes())
    s = _CACHE.get(key)
    if s is None:
        s = SOEStructure(G, int(ham.n_data))
        _CACHE[key] = s
    return s


def _physical(S: SOEStructure, ham, u: np.ndarray) -> np.ndarray:
    """e_L + u G for a batch u [R, m] (uint8 in, uint8 out), via an exact float64 matmul."""
    uf = np.atleast_2d(np.asarray(u)).astype(np.float64)
    e = (np.rint(uf @ S.G_f).astype(np.int64) & 1).astype(np.uint8)
    return e ^ ham.seed[None, :]


def _delta(ex: np.ndarray, ez: np.ndarray, q: np.ndarray, fx: np.ndarray, fz: np.ndarray) -> np.ndarray:
    exq, ezq = ex[q], ez[q]
    old = (exq | ezq).astype(np.int16)
    new = ((exq ^ fx) | (ezq ^ fz)).astype(np.int16)
    return (new - old).sum(axis=1).astype(np.int64)


def class_histogram(pool_u: np.ndarray, ham, radius: int):
    """Return (E0, N) with N[t] = number of ball states of Pauli weight E0 + t - OFF."""
    S = get_structure(ham)
    n = S.n
    phys = _physical(S, ham, pool_u)
    wt = np.count_nonzero(phys[:, :n] | phys[:, n:], axis=1)
    c0 = int(np.argsort(wt, kind="stable")[0])          # same centre as adaptive_soe._score_one_class
    E0 = int(wt[c0])
    ex, ez = phys[c0, :n], phys[c0, n:]
    d1 = _delta(ex, ez, S.q1, S.fx1, S.fz1)              # [m]
    N = np.zeros(S.L, dtype=np.int64)
    N[S.OFF] += 1
    N += np.bincount(d1 + S.OFF, minlength=S.L)
    if radius >= 2:
        K = S.K
        h = np.bincount(d1 + K, minlength=2 * K + 1)      # index i <-> dE = i - K
        conv = np.convolve(h, h)                          # index s <-> dE_j + dE_k = s - 2K  (ordered, incl. j == k)
        diag = np.zeros_like(conv)
        diag[2 * np.arange(2 * K + 1)] = h                # j == k terms
        N += (conv - diag) // 2                           # unordered j < k, all treated as disjoint
        if S.n_overlap_pairs:
            naive = d1[S.pj] + d1[S.pk]
            true = _delta(ex, ez, S.q2, S.fx2, S.fz2)
            N -= np.bincount(naive + S.OFF, minlength=S.L)
            N += np.bincount(true + S.OFF, minlength=S.L)
    return E0, N


def _expected_count(m: int, radius: int) -> int:
    return 1 + m + (m * (m - 1) // 2 if radius >= 2 else 0)


def fast_soe_prediction(pools, hams, p, top_m=1, radius=1, return_hist=False):
    """Drop-in replacement for adaptive_soe.fixed_soe_prediction(pools, hams, p, 1, radius), radius in {1, 2}.

    Returns (pred, scores, counts, minima) like the original; with return_hist=True also the
    per-class (E0, N) histograms (used by the consistency check to identify exact ties).
    """
    if int(top_m) != 1 or int(radius) not in (1, 2):
        raise ValueError("fast_soe_prediction supports top_m=1 and radius 1 or 2 only")
    w = (float(p) / 3.0) / (1.0 - float(p))
    scores = np.zeros(4, dtype=float)
    counts = np.zeros(4, dtype=int)
    minima = np.zeros(4, dtype=int)
    hists = {}
    for c in range(4):
        E0, N = class_histogram(pools[c], hams[c], int(radius))
        S = get_structure(hams[c])
        if N.min() < 0 or int(N.sum()) != _expected_count(S.m, int(radius)):
            raise RuntimeError(f"SOE histogram inconsistent (class {c}): min {N.min()}, sum {N.sum()}")
        nz = np.flatnonzero(N)
        expo = E0 + nz - S.OFF                            # ascending Pauli weights
        if expo.min() < 0:
            raise RuntimeError("negative Pauli weight in SOE histogram")
        scores[c] = float(np.sum(N[nz].astype(np.float64) * w ** expo.astype(np.float64)))
        counts[c] = int(N.sum())
        minima[c] = int(expo.min())
        if return_hist:
            hists[c] = (tuple(expo.tolist()), tuple(N[nz].tolist()))
    pred = _stable_argmax(scores)
    if return_hist:
        return pred, scores, counts, minima, hists
    return pred, scores, counts, minima
