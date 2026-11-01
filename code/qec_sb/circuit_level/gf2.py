"""GF(2) linear algebra on Python-int bitsets (rows as ints, bit i = column i).

All matrices are "list of int rows" with an explicit column count. This is fast
enough for DEM sizes (thousands of detectors x tens of thousands of mechanisms).
"""
from __future__ import annotations

import numpy as np


# ------------------------------------------------------------------ conversions
def rows_from_dense(M: np.ndarray) -> list[int]:
    """Dense 0/1 array (r, c) -> list of int rows."""
    M = np.asarray(M, dtype=np.uint8)
    out = []
    for row in M:
        packed = np.packbits(row, bitorder="little").tobytes()
        out.append(int.from_bytes(packed, "little"))
    return out


def dense_from_rows(rows: list[int], ncols: int) -> np.ndarray:
    nbytes = (ncols + 7) // 8
    out = np.zeros((len(rows), ncols), dtype=np.uint8)
    for k, r in enumerate(rows):
        b = np.frombuffer(r.to_bytes(nbytes, "little"), dtype=np.uint8)
        out[k] = np.unpackbits(b, bitorder="little")[:ncols]
    return out


def rows_from_csr(H) -> list[int]:
    """scipy CSR (0/1) -> int rows."""
    H = H.tocsr()
    out = []
    for i in range(H.shape[0]):
        v = 0
        for j in H.indices[H.indptr[i]:H.indptr[i + 1]]:
            v ^= 1 << int(j)
        out.append(v)
    return out


def popcount(x: int) -> int:
    return x.bit_count()          # Python >= 3.10


def support(x: int) -> list[int]:
    out = []
    while x:
        low = x & -x
        out.append(low.bit_length() - 1)
        x ^= low
    return out


# ------------------------------------------------------------------ elimination
class Basis:
    """Incremental GF(2) basis with pivot = highest set bit. Supports reduce/insert."""

    def __init__(self):
        self.piv: dict[int, int] = {}

    def reduce(self, v: int) -> int:
        piv = self.piv
        while v:
            c = v.bit_length() - 1
            p = piv.get(c)
            if p is None:
                return v
            v ^= p
        return 0

    def insert(self, v: int) -> bool:
        """Insert v if independent. Returns True if rank increased."""
        v = self.reduce(v)
        if v == 0:
            return False
        self.piv[v.bit_length() - 1] = v
        return True

    @property
    def rank(self) -> int:
        return len(self.piv)


def rank(rows: list[int]) -> int:
    b = Basis()
    for r in rows:
        b.insert(r)
    return b.rank


def rref(rows: list[int]) -> dict[int, int]:
    """Fully reduced row echelon form: pivot_col -> row (no other pivot bits set)."""
    b = Basis()
    for r in rows:
        b.insert(r)
    piv = b.piv
    for c in sorted(piv):                      # ascending: rows with higher pivots may hold bit c
        pr = piv[c]
        for c2 in piv:
            if c2 > c and (piv[c2] >> c) & 1:
                piv[c2] ^= pr
    return piv


def nullspace(rows: list[int], ncols: int) -> list[int]:
    """Basis of {x : rows . x = 0}. Returns list of int vectors over ncols bits."""
    piv = rref(rows)
    if not piv:
        return [1 << f for f in range(ncols)]
    pcols = sorted(piv)
    pset = set(pcols)
    free = [c for c in range(ncols) if c not in pset]
    # R: (rank, ncols) dense; kernel vector for free f: e_f + sum_{pivot c: R[c,f]=1} e_c
    R = dense_from_rows([piv[c] for c in pcols], ncols)          # rows ordered as pcols
    Rf = R[:, free]                                               # (rank, nfree)
    pc = np.asarray(pcols)
    kernel = []
    for k, f in enumerate(free):
        v = 1 << f
        for c in pc[Rf[:, k].astype(bool)]:
            v |= 1 << int(c)
        kernel.append(v)
    return kernel


def solve(rows: list[int], rhs: np.ndarray, ncols: int) -> int | None:
    """One solution x of rows . x = rhs (GF(2)), or None. rhs: 0/1 array len(rows)."""
    # augment each row with its rhs bit at position 0 (the LOWEST bit): pivots are highest set bits,
    # so the rhs column must never become a pivot. Shift the matrix columns up by one.
    aug = [(r << 1) | int(b) for r, b in zip(rows, rhs)]
    piv = rref(aug)
    if 0 in piv:                           # a row reduced to [0 ... 0 | 1]: inconsistent
        return None
    x = 0
    for c, r in piv.items():               # fully reduced: x_c = rhs bit of the pivot row (free vars = 0)
        if r & 1:
            x |= 1 << (c - 1)
    return x


def complement_basis(sub_rows: list[int], big_rows: list[int]) -> list[int]:
    """Vectors of big_rows that extend span(sub_rows); returned reduced against the sub basis."""
    b = Basis()
    for r in sub_rows:
        b.insert(r)
    out = []
    for r in big_rows:
        v = b.reduce(r)
        if v:
            b.piv[v.bit_length() - 1] = v
            out.append(v)
    return out


# ------------------------------------------------------------------ low-weight search
def low_weight_vectors(gens: list[int], ncols: int, trials: int, rng: np.random.Generator,
                       max_weight: int | None = None) -> dict[int, int]:
    """Randomised information-set search for low-weight vectors in span(gens).

    Repeats: random column order -> echelon form with pivots chosen in that order
    -> every resulting row is a span element; keep the light ones. Returns
    {vector: weight}. Standard trick for finding sparse generators of an LDPC-like space.
    """
    found: dict[int, int] = {}
    if not gens:
        return found
    for _ in range(trials):
        perm = rng.permutation(ncols)
        # relabel bits by perm so that "highest bit" pivoting follows a random order
        pos = np.empty(ncols, dtype=np.int64)
        pos[perm] = np.arange(ncols)
        remapped = [_remap(g, pos, ncols) for g in gens]
        piv = rref(remapped)
        inv = perm            # bit k in remapped space = column perm[k]
        for r in piv.values():
            w = popcount(r)
            if max_weight is not None and w > max_weight:
                continue
            v = _remap(r, inv, ncols)
            if v not in found:
                found[v] = w
    return found


def _remap(v: int, table: np.ndarray, ncols: int) -> int:
    out = 0
    for i in support(v):
        out |= 1 << int(table[i])
    return out
