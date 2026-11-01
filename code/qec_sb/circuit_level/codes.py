"""Code constructions: bivariate-bicycle (BB) qLDPC codes and their logical operators."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import gf2


@dataclass(frozen=True)
class CSSCode:
    name: str
    hx: np.ndarray          # (mx, n) uint8
    hz: np.ndarray          # (mz, n) uint8
    lx: np.ndarray          # (k, n)  X logicals
    lz: np.ndarray          # (k, n)  Z logicals

    @property
    def n(self):
        return self.hx.shape[1]

    @property
    def k(self):
        return self.lz.shape[0]


# ---------------------------------------------------------------- BB codes
BB_FAMILY = {
    # name: (l, m, A terms, B terms); terms as (xpow, ypow). Bravyi et al., Nature 627, 778 (2024)
    "bb72":  (6, 6,  [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)]),     # [[72,12,6]]
    "bb90":  (15, 3, [(9, 0), (0, 1), (0, 2)], [(0, 0), (2, 0), (7, 0)]),     # [[90,8,10]]
    "bb108": (9, 6,  [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)]),     # [[108,8,10]]
    "bb144": (12, 6, [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)]),     # [[144,12,12]]
}


def _shift(size, power):
    return np.roll(np.eye(size, dtype=np.uint8), power, axis=1)


def _poly(l, m, terms):
    M = np.zeros((l * m, l * m), dtype=np.uint8)
    for a, b in terms:
        M ^= np.kron(_shift(l, a), _shift(m, b))
    return M


def bb_monomials(name: str):
    """The individual monomial permutation matrices: A = A1 + A2 + A3, B = B1 + B2 + B3 (Bravyi et al. ordering:
    for bb72/bb144 A1 = x^3, A2 = y, A3 = y^2, B1 = y^3, B2 = x, B3 = x^2). M[r, c] = 1 iff c = r + monomial."""
    l, m, aterms, bterms = BB_FAMILY[name]
    A = [np.kron(_shift(l, a), _shift(m, b)) for a, b in aterms]
    B = [np.kron(_shift(l, a), _shift(m, b)) for a, b in bterms]
    return A, B


def bb_code(name: str) -> CSSCode:
    l, m, aterms, bterms = BB_FAMILY[name]
    A, B = _poly(l, m, aterms), _poly(l, m, bterms)
    hx = np.concatenate([A, B], axis=1)
    hz = np.concatenate([B.T, A.T], axis=1)
    assert not ((hx.astype(int) @ hz.T.astype(int)) % 2).any(), "HX HZ^T != 0"
    lx, lz = css_logicals(hx, hz)
    return CSSCode(name, hx, hz, lx, lz)


# ---------------------------------------------------------------- logicals
def css_logicals(hx: np.ndarray, hz: np.ndarray):
    """Z logicals: ker(hx) modulo rowspace(hz). X logicals: ker(hz) modulo rowspace(hx).

    Returned pairs are made symplectically conjugate (lx[i] . lz[j] = delta_ij).
    """
    n = hx.shape[1]
    kz = gf2.nullspace(gf2.rows_from_dense(hx), n)             # candidates for Z logicals
    lz_rows = gf2.complement_basis(gf2.rows_from_dense(hz), kz)
    kx = gf2.nullspace(gf2.rows_from_dense(hz), n)
    lx_rows = gf2.complement_basis(gf2.rows_from_dense(hx), kx)
    assert len(lz_rows) == len(lx_rows), (len(lz_rows), len(lx_rows))
    lz = gf2.dense_from_rows(lz_rows, n)
    lx = gf2.dense_from_rows(lx_rows, n)
    lx = _symplectic_pair(lx, lz)
    return lx, lz


def _symplectic_pair(lx, lz):
    """Row-reduce lx so that lx @ lz.T = I (mod 2)."""
    lx = lx.copy()
    k = lz.shape[0]
    for i in range(k):
        ov = (lx.astype(int) @ lz[i].astype(int)) % 2
        js = [j for j in range(i, k) if ov[j]]
        assert js, "logical pairing failed (lz not independent of lx?)"
        j = js[0]
        if j != i:
            lx[[i, j]] = lx[[j, i]]
            ov[[i, j]] = ov[[j, i]]
        for j2 in range(k):
            if j2 != i and ov[j2]:
                lx[j2] ^= lx[i]
    return lx


def min_weight_estimate(rows_dense: np.ndarray, mod_rows: np.ndarray, trials=200, seed=0):
    """Cheap estimate of min weight of a coset space (e.g. logicals modulo stabilizers)."""
    n = rows_dense.shape[1]
    rng = np.random.default_rng(seed)
    gens = gf2.rows_from_dense(np.concatenate([rows_dense, mod_rows], axis=0))
    stab = gf2.Basis()
    for r in gf2.rows_from_dense(mod_rows):
        stab.insert(r)
    found = gf2.low_weight_vectors(gens, n, trials, rng)
    best = None
    for v, w in found.items():
        if stab.reduce(v) and (best is None or w < best):
            best = w
    return best
