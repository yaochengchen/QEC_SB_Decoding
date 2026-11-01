from dataclasses import dataclass
from itertools import combinations
import os
import numpy as np

from .rotated_code import RotatedSurfaceCode

# How single-generator flip energies dE_j are evaluated (used by CP-SB here and by CP-SA in cp_sa.py):
#   "sparse" (default): only the <= K qubits that generator j touches are read  -> O(m K) per replica and step
#   "dense"           : the original evaluation on the full 2n-bit error vector  -> O(m n)
# Both give the same integers, so trajectories, pools and decisions are identical (checked by
# scripts/run_soe_fast.py check). Override with the environment variable QECSB_DELTA=dense.
DELTA_MODE = os.environ.get("QECSB_DELTA", "sparse")
_SUPPORT_CACHE: dict = {}


class GeneratorSupport:
    """Per-code padded supports of the generators: qubit index and X/Z flip bit on each touched qubit."""

    def __init__(self, G: np.ndarray, n: int):
        G = np.asarray(G, dtype=np.uint8)
        m = int(G.shape[0])
        gx, gz = G[:, :n], G[:, n:]
        supp = (gx | gz).astype(bool)
        K = max(1, int(supp.sum(axis=1).max())) if m else 1
        self.m, self.n, self.K = m, int(n), K
        self.G_f = G.astype(np.float64)
        # padding: qubit 0 with zero flips, which contributes exactly 0 to every dE
        self.q = np.zeros((m, K), dtype=np.int64)
        self.fx = np.zeros((m, K), dtype=np.uint8)
        self.fz = np.zeros((m, K), dtype=np.uint8)
        self.lists = []
        for j in range(m):
            qs = np.flatnonzero(supp[j])
            self.q[j, :len(qs)] = qs
            self.fx[j, :len(qs)] = gx[j, qs]
            self.fz[j, :len(qs)] = gz[j, qs]
            self.lists.append(tuple((int(q), int(gx[j, q]), int(gz[j, q])) for q in qs))


def generator_support(ham) -> GeneratorSupport:
    G = ham.generators
    key = (int(ham.n_data), G.shape, G.tobytes())
    s = _SUPPORT_CACHE.get(key)
    if s is None:
        s = GeneratorSupport(G, int(ham.n_data))
        _SUPPORT_CACHE[key] = s
    return s


def sparse_flip_deltas(sup: GeneratorSupport, seed: np.ndarray, ubits: np.ndarray) -> np.ndarray:
    """dE_j for every generator j and every row of ubits [R, m] (int64 [R, m]), reading only supp(j)."""
    n = sup.n
    phys = (np.rint(ubits.astype(np.float64) @ sup.G_f).astype(np.int64) & 1).astype(np.uint8) ^ seed[None, :]
    exq, ezq = phys[:, :n][:, sup.q], phys[:, n:][:, sup.q]          # [R, m, K]
    new = ((exq ^ sup.fx) | (ezq ^ sup.fz)).astype(np.int16)
    old = (exq | ezq).astype(np.int16)
    return (new - old).sum(axis=2).astype(np.int64)


@dataclass(frozen=True)
class StabilizerHamiltonian:
    """Exact Pauli-weight Hamiltonian in stabilizer-generator spin variables.

    For a fixed syndrome and logical class choose a feasible seed error e0. Every
    error in that class is e(u)=e0 XOR (u @ G), where rows of G are physical
    stabilizer masks and u are binary generator-selection variables.

    With spins s_j=(-1)^u_j, the depolarizing Pauli weight is represented exactly
    at spin vertices by a multilinear polynomial
        E(s) = sum_A c_A prod_{j in A} s_j.
    On a rotated surface code a data qubit touches at most two X and two Z checks,
    so terms are at most four-body. There is no syndrome penalty and no slack bit.
    """

    generators: np.ndarray      # [m, 2n] GF(2)
    seed: np.ndarray            # [2n] GF(2)
    terms: tuple                # tuple[(indices_tuple, coefficient), ...]
    n_data: int

    @property
    def n_spins(self):
        return int(self.generators.shape[0])

    def physical_from_bits(self, u):
        u = np.asarray(u, dtype=np.uint8)
        one = u.ndim == 1
        uu = u[None, :] if one else u
        masks = (uu @ self.generators) & 1
        out = masks ^ self.seed[None, :]
        return out[0] if one else out

    def physical_from_spins(self, s):
        s = np.asarray(s)
        u = (s < 0).astype(np.uint8)
        return self.physical_from_bits(u)

    def energy_bits(self, u):
        """Exact integer Pauli weight for binary states [m] or batches [R,m]."""
        u = np.asarray(u)
        if u.ndim not in (1, 2) or u.shape[-1] != self.n_spins:
            raise ValueError("u must have shape [n_spins] or [batch, n_spins]")
        if not np.all((u == 0) | (u == 1)):
            raise ValueError("u must contain only binary values 0 or 1")
        physical = self.physical_from_bits(u)
        one = physical.ndim == 1
        batch = physical[None, :] if one else physical
        weights = np.count_nonzero(
            batch[:, :self.n_data] | batch[:, self.n_data:], axis=1
        )
        return int(weights[0]) if one else weights.astype(int)

    def energy(self, s):
        """Multilinear energy. Equals physical Pauli weight for s in {+/-1}."""
        s = np.asarray(s, dtype=float)
        one = s.ndim == 1
        ss = s[None, :] if one else s
        e = np.zeros(ss.shape[0], dtype=float)
        for inds, coeff in self.terms:
            if len(inds) == 0:
                e += coeff
            else:
                e += coeff * np.prod(ss[:, inds], axis=1)
        return e[0] if one else e

    def minimizing_force(self, source):
        """Return -grad E for a multilinear source array [R,m]."""
        source = np.asarray(source, dtype=float)
        one = source.ndim == 1
        x = source[None, :] if one else source
        force = np.zeros_like(x)
        for inds, coeff in self.terms:
            if not inds:
                continue
            inds = tuple(inds)
            if len(inds) == 1:
                force[:, inds[0]] -= coeff
                continue
            for j in inds:
                others = [k for k in inds if k != j]
                force[:, j] -= coeff * np.prod(x[:, others], axis=1)
        return force[0] if one else force

    def force_l1_bound(self):
        b = np.zeros(self.n_spins, dtype=float)
        for inds, coeff in self.terms:
            for j in inds:
                b[j] += abs(coeff)
        return b


def stabilizer_generator_matrix(code: RotatedSurfaceCode):
    """Physical (X|Z) masks for an independent CSS stabilizer generating set."""
    n = code.n
    z0 = np.zeros(n, dtype=np.uint8)
    rows = []
    # X-type stabilizers toggle the X component of an error.
    for row in code.Hx:
        rows.append(np.concatenate([row.astype(np.uint8), z0]))
    # Z-type stabilizers toggle the Z component.
    for row in code.Hz:
        rows.append(np.concatenate([z0, row.astype(np.uint8)]))
    return np.asarray(rows, dtype=np.uint8)


def build_stabilizer_hamiltonian(code: RotatedSurfaceCode, seed_ex, seed_ez):
    """Build exact stabilizer-variable Pauli-weight polynomial for one logical class."""
    n = code.n
    seed = np.concatenate([
        np.asarray(seed_ex, dtype=np.uint8),
        np.asarray(seed_ez, dtype=np.uint8),
    ])
    G = stabilizer_generator_matrix(code)

    # Aggregate identical monomials. For one qubit:
    #   1[Pauli != I] = 3/4 - tx/4 - tz/4 - tx*tz/4,
    # where tx=(-1)^ex and tz=(-1)^ez.
    coeff = {(): 0.75 * n}
    for q in range(n):
        sx = tuple(np.flatnonzero(G[:, q]).tolist())
        sz = tuple(np.flatnonzero(G[:, n + q]).tolist())
        sign_x = 1.0 if seed[q] == 0 else -1.0
        sign_z = 1.0 if seed[n + q] == 0 else -1.0
        cross = tuple(sorted(set(sx) ^ set(sz)))  # disjoint for CSS, XOR is robust.
        for inds, c in (
            (sx, -0.25 * sign_x),
            (sz, -0.25 * sign_z),
            (cross, -0.25 * sign_x * sign_z),
        ):
            coeff[inds] = coeff.get(inds, 0.0) + c

    terms = tuple((inds, float(c)) for inds, c in sorted(coeff.items()) if abs(c) > 1e-15)
    return StabilizerHamiltonian(G, seed, terms, n)


@dataclass
class CPSBParams:
    replicas: int = 64
    steps: int = 500
    dt: float = 0.45
    coupling_scale: float = 0.70
    warm_strength: float = 0.08
    random_sigma: float = 0.01
    individualized: bool = False
    gb_control_A: float = 0.20


def run_constraint_preserving_sb(ham: StabilizerHamiltonian, params: CPSBParams, seed: int,
                                 warm_bits=None, return_trace=False):
    """Discrete simulated-bifurcation dynamics on stabilizer variables.

    The force is the exact multilinear -grad(E) evaluated on sign(x), analogous
    to dSB's discrete interaction field. Any returned generator bitstring maps to
    a syndrome-feasible correction in the selected logical class by construction.

    If individualized=True, an individualized restoring coefficient inspired by
    generalized SB is used. This is optional; the core CP property is independent
    of it. When ``return_trace`` is true, return ``(spins, trace)`` with physical
    Pauli-weight summaries recorded after every iteration. Trace collection does
    not draw randomness or modify the dynamics.
    """
    rng = np.random.default_rng(seed)
    R, m = int(params.replicas), ham.n_spins
    x = rng.uniform(-0.1, 0.1, size=(R, m))
    y = np.zeros_like(x)
    if warm_bits is not None:
        wb = np.asarray(warm_bits, dtype=np.uint8)
        warm_spin = 1.0 - 2.0 * wb.astype(float)
        x = params.warm_strength * warm_spin[None, :] + rng.normal(
            0.0, params.random_sigma, size=(R, m)
        )

    # Normalize the nonlinear force by a deterministic per-variable L1 bound.
    # This keeps dt/coupling settings comparable across syndromes/classes.
    bound = ham.force_l1_bound()
    norm = max(float(np.max(bound)), 1.0)
    c = float(params.coupling_scale) / norm
    p_ind = np.ones((R, m), dtype=float)

    # d=3 has only 8 stabilizer variables. Precompute the exact discrete force
    # for every sign vertex to make repeated SB time steps cheap. For larger codes
    # fall back to direct multilinear-force evaluation.
    force_table = None
    powers = None
    if m <= 16:
        ids = np.arange(1 << m, dtype=np.uint32)[:, None]
        bits = ((ids >> np.arange(m, dtype=np.uint32)) & 1).astype(np.uint8)
        vertices = 1.0 - 2.0 * bits.astype(float)
        force_table = ham.minimizing_force(vertices)
        powers = (1 << np.arange(m, dtype=np.uint32)).astype(np.uint32)

    sparse = force_table is None and DELTA_MODE != "dense"
    sup = generator_support(ham) if sparse else None

    M = int(params.steps)
    dt = float(params.dt)
    trace = [] if return_trace else None
    for step in range(M):
        discrete = np.sign(x)
        discrete[discrete == 0] = 1.0
        if force_table is not None:
            state_ids = ((discrete < 0).astype(np.uint32) @ powers).astype(np.int64)
            force = force_table[state_ids]
        elif sparse:
            # Same integers dE_j as the dense branch below, reading only the qubits of generator j.
            delta = sparse_flip_deltas(sup, ham.seed, (discrete < 0).astype(np.uint8)).astype(float)
            force = 0.5 * delta * discrete
        else:
            # Exact discrete multilinear force via physical Pauli-weight changes.
            # Flipping stabilizer spin j toggles physical mask G_j. At a sign
            # vertex, DeltaE_j = E(s_j->-s_j)-E(s) = -2*s_j*dE/ds_j, hence
            # the minimizing force -dE/ds_j = 0.5*DeltaE_j*s_j. This avoids
            # Python loops over monomials and is much faster for d>=5.
            ubits = (discrete < 0).astype(np.uint8)
            phys = ham.physical_from_bits(ubits)
            n = ham.n_data
            ex = phys[:, :n]
            ez = phys[:, n:]
            gx = ham.generators[:, :n]
            gz = ham.generators[:, n:]
            old_nz = (ex | ez)[:, None, :]
            new_nz = ((ex[:, None, :] ^ gx[None, :, :]) |
                      (ez[:, None, :] ^ gz[None, :, :]))
            delta = new_nz.sum(axis=2).astype(float) - old_nz.sum(axis=2).astype(float)
            force = 0.5 * delta * discrete
        if params.individualized:
            denom = max(M - step, 1)
            p_ind -= (1.0 - params.gb_control_A * x * x) * p_ind / denom
            restoring = -p_ind * x
        else:
            a = (step + 1) / M
            restoring = -(1.0 - a) * x
        y += (restoring + c * force) * dt
        x += y * dt
        hit = np.abs(x) > 1.0
        x = np.where(hit, np.sign(x), x)
        y = np.where(hit, 0.0, y)

        if return_trace:
            trace_spins = np.sign(x)
            trace_spins[trace_spins == 0] = 1.0
            physical = ham.physical_from_spins(trace_spins)
            weights = np.count_nonzero(
                physical[:, :ham.n_data] | physical[:, ham.n_data:], axis=1
            )
            trace.append({
                "iteration": step + 1,
                "best_weight": float(np.min(weights)),
                "mean_weight": float(np.mean(weights)),
                "feasible_fraction": 1.0,
            })

    s = np.sign(x)
    s[s == 0] = 1.0
    spins = s.astype(np.int8)
    return (spins, trace) if return_trace else spins


def stabilizer_hamming_neighborhood(center_bits, n_spins, radius):
    """Enumerate generator bitstrings within Hamming radius around center."""
    center = np.asarray(center_bits, dtype=np.uint8)
    out = [center.copy()]
    for r in range(1, int(radius) + 1):
        for comb in combinations(range(n_spins), r):
            v = center.copy()
            v[list(comb)] ^= 1
            out.append(v)
    return np.asarray(out, dtype=np.uint8)
