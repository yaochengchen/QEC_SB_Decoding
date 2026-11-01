from dataclasses import dataclass
import math
import numpy as np
from .rotated_code import RotatedSurfaceCode

@dataclass
class QuboModel:
    S: np.ndarray              # symmetric matrix: E(v)=v^T S v + const
    slack_specs: list
    physical_n: int
    pauli_weight: float
    penalty: float
    constant: float


def _add_square(S, coeff, rhs, lam):
    nz = np.flatnonzero(coeff)
    for i in nz:
        ai = coeff[i]
        S[i, i] += lam * (ai * ai - 2.0 * rhs * ai)
    for pos, i in enumerate(nz):
        ai = coeff[i]
        for j in nz[pos + 1:]:
            value = lam * ai * coeff[j]
            S[i, j] += value
            S[j, i] += value


def build_depolarizing_qubo(code: RotatedSurfaceCode, sx, sz, p: float, penalty_mult=None) -> QuboModel:
    """QUBO for correlated depolarizing MAP-configuration decoding.

    Physical bits are x_i,z_i. For depolarizing noise I/X/Z/Y have costs 0/w/w/w,
    represented exactly by w * (x_i + z_i - x_i z_i).

    Syndrome parity constraints are quadratized using integer slack variables:
        sum(bits on check) - 2 k = syndrome_bit.
    """
    n = code.n
    constraints = []
    for row, b in zip(code.Hz, sx):
        constraints.append((np.flatnonzero(row).tolist(), int(b), "x"))
    for row, b in zip(code.Hx, sz):
        constraints.append((np.flatnonzero(row).tolist(), int(b), "z"))

    offset = 2 * n
    slack_specs = []
    for phys, _, _ in constraints:
        max_k = len(phys) // 2
        nbits = max(1, math.ceil(math.log2(max_k + 1))) if max_k else 0
        weights = [1 << j for j in range(nbits)]
        inds = list(range(offset, offset + nbits))
        offset += nbits
        slack_specs.append((inds, weights))

    S = np.zeros((offset, offset), dtype=float)
    w = math.log((1.0 - p) / (p / 3.0))
    for i in range(n):
        S[i, i] += w
        S[n + i, n + i] += w
        # symmetric matrix => pair coefficient is 2*S[i,j]
        S[i, n + i] -= w / 2.0
        S[n + i, i] -= w / 2.0

    if penalty_mult is None:
        # Conservative generic setting. For d=3 experiments we certify and use 5.
        penalty_mult = n + 1
    lam = float(penalty_mult) * w

    constant = 0.0
    for (phys, rhs, kind), (inds, weights) in zip(constraints, slack_specs):
        coeff = np.zeros(offset, dtype=float)
        base = 0 if kind == "x" else n
        for q in phys:
            coeff[base + q] = 1.0
        for si, sw in zip(inds, weights):
            coeff[si] = -2.0 * sw
        _add_square(S, coeff, rhs, lam)
        constant += lam * rhs * rhs

    return QuboModel(S=S, slack_specs=slack_specs, physical_n=n, pauli_weight=w, penalty=lam, constant=constant)


def qubo_to_zero_field_ising(S: np.ndarray):
    """Convert symmetric QUBO to zero-field Ising using one gauge/auxiliary spin.

    Binary v=(1-s)/2. Local fields are coupled to an auxiliary spin. At readout,
    physical spins are interpreted relative to the auxiliary spin.
    """
    N = S.shape[0]
    h = 0.5 * (S @ np.ones(N))
    J = -0.5 * S.copy()
    np.fill_diagonal(J, 0.0)
    aug = np.zeros((N + 1, N + 1), dtype=float)
    aug[:N, :N] = J
    aug[:N, N] = h
    aug[N, :N] = h
    scale = max(float(np.max(np.abs(aug))), 1e-12)
    return aug / scale, scale


def warm_binary_from_correction(code, model: QuboModel, ex, ez, sx, sz):
    n = code.n
    v = np.zeros(model.S.shape[0], dtype=np.uint8)
    v[:n] = ex
    v[n:2*n] = ez
    constraints = []
    for row, b in zip(code.Hz, sx):
        constraints.append((np.flatnonzero(row), int(b), "x"))
    for row, b in zip(code.Hx, sz):
        constraints.append((np.flatnonzero(row), int(b), "z"))
    for (phys, rhs, kind), (inds, weights) in zip(constraints, model.slack_specs):
        arr = ex if kind == "x" else ez
        k = (int(arr[phys].sum()) - rhs) // 2
        if k < 0:
            raise ValueError("Warm correction is not syndrome-feasible.")
        for ind, weight in zip(inds, weights):
            v[ind] = (k // weight) & 1
    return v
