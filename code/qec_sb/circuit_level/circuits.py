"""Stim circuits.

* rotated surface code memory: ``stim.Circuit.generated`` (standard).
* BB codes: a generic circuit-level memory-Z experiment with a greedy edge-colouring
  CNOT schedule and uniform circuit-level depolarising noise.

  (generic schedule, NOT used for paper 1; the BB row uses bb_memory_z_bravyi below)
"""
from __future__ import annotations

import numpy as np

from .codes import CSSCode


def surface_memory(d: int, p: float, rounds: int | None = None, basis="z"):
    import stim
    rounds = rounds or d
    return stim.Circuit.generated(
        f"surface_code:rotated_memory_{basis}", distance=d, rounds=rounds,
        after_clifford_depolarization=p, before_round_data_depolarization=p,
        before_measure_flip_probability=p, after_reset_flip_probability=p)


# ---------------------------------------------------------------- edge colouring
def edge_colouring(h: np.ndarray) -> list[list[tuple[int, int]]]:
    """Greedy proper edge colouring of the bipartite Tanner graph of h (checks x qubits).

    Returns layers; each layer is a list of (check, qubit) edges with no shared endpoint.
    Greedy with Vizing-style repair is overkill here: for degree-6 BB graphs greedy
    typically needs 7-8 colours; depth is reported by the circuit builder.
    """
    checks, qubits = np.nonzero(h)
    edges = list(zip(checks.tolist(), qubits.tolist()))
    # order edges by (check, qubit) degree to help greedy
    deg_c = h.sum(1).astype(int)
    deg_q = h.sum(0).astype(int)
    edges.sort(key=lambda e: -(deg_c[e[0]] + deg_q[e[1]]))
    colour_c: dict[int, set] = {}
    colour_q: dict[int, set] = {}
    layers: list[list[tuple[int, int]]] = []
    for c, q in edges:
        used = colour_c.setdefault(c, set()) | colour_q.setdefault(q, set())
        k = 0
        while k in used:
            k += 1
        while len(layers) <= k:
            layers.append([])
        layers[k].append((c, q))
        colour_c[c].add(k)
        colour_q[q].add(k)
    return layers


# ---------------------------------------------------------------- BB memory circuit
def bb_memory_z(code: CSSCode, p: float, rounds: int):
    """Memory-Z experiment: data |0>, `rounds` rounds of X and Z syndrome extraction, final Z readout.

    Qubit layout: data 0..n-1, X-checks n..n+mx-1, Z-checks n+mx..n+mx+mz-1.
    Noise (uniform p): DEPOLARIZE2 after every CNOT, DEPOLARIZE1 on idle qubits in every CNOT
    layer, X/Z flip after reset, flip before measurement, DEPOLARIZE1 on data before each round.
    Detectors: Z-checks every round (first round vs |0> init), X-checks from round 2 on,
    final data readout vs last Z-check round. Observables: the k logical Z operators.
    """
    import stim
    n, mx, mz = code.n, code.hx.shape[0], code.hz.shape[0]
    data = list(range(n))
    xq = [n + i for i in range(mx)]
    zq = [n + mx + i for i in range(mz)]
    allq = data + xq + zq
    x_layers = edge_colouring(code.hx)
    z_layers = edge_colouring(code.hz)

    c = stim.Circuit()
    for q in allq:
        c.append("QUBIT_COORDS", [q], [q])
    meas_index = {}        # (kind, check, round) -> absolute measurement index
    counter = [0]

    def measure(kind, qubits, keys, mbasis):
        c.append("X_ERROR" if mbasis == "z" else "Z_ERROR", qubits, p)
        c.append("M" if mbasis == "z" else "MX", qubits)
        for k in keys:
            meas_index[k] = counter[0]
            counter[0] += 1

    def rec(key):
        return stim.target_rec(meas_index[key] - counter[0])

    def cnot_layer(pairs):
        flat = [q for pair in pairs for q in pair]
        c.append("CX", flat)
        c.append("DEPOLARIZE2", flat, p)
        busy = set(flat)
        idle = [q for q in allq if q not in busy]
        if idle:
            c.append("DEPOLARIZE1", idle, p)
        c.append("TICK")

    c.append("R", data)
    c.append("X_ERROR", data, p)
    c.append("TICK")
    for r in range(rounds):
        c.append("DEPOLARIZE1", data, p)
        c.append("RX", xq)
        c.append("Z_ERROR", xq, p)
        c.append("R", zq)
        c.append("X_ERROR", zq, p)
        c.append("TICK")
        for layer in x_layers:                          # X-check -> data
            cnot_layer([(xq[chk], q) for chk, q in layer])
        for layer in z_layers:                          # data -> Z-check
            cnot_layer([(q, zq[chk]) for chk, q in layer])
        measure("x", xq, [("x", i, r) for i in range(mx)], "x")
        measure("z", zq, [("z", i, r) for i in range(mz)], "z")
        for i in range(mz):
            if r == 0:
                c.append("DETECTOR", [rec(("z", i, r))], [i, r, 0])
            else:
                c.append("DETECTOR", [rec(("z", i, r)), rec(("z", i, r - 1))], [i, r, 0])
        if r > 0:
            for i in range(mx):
                c.append("DETECTOR", [rec(("x", i, r)), rec(("x", i, r - 1))], [i, r, 1])
        c.append("TICK")
    measure("d", data, [("d", q, 0) for q in data], "z")
    for i in range(mz):
        targets = [rec(("d", q, 0)) for q in np.nonzero(code.hz[i])[0]] + [rec(("z", i, rounds - 1))]
        c.append("DETECTOR", targets, [i, rounds, 0])
    for a in range(code.k):
        targets = [rec(("d", q, 0)) for q in np.nonzero(code.lz[a])[0]]
        c.append("OBSERVABLE_INCLUDE", targets, a)
    return c


def circuit_depth_report(code: CSSCode) -> dict:
    return {"x_cnot_layers": len(edge_colouring(code.hx)), "z_cnot_layers": len(edge_colouring(code.hz))}


# ---------------------------------------------------------------- Bravyi et al. depth-8 BB circuit
# Syndrome cycle of Bravyi et al., Nature 627, 778 (2024) (github.com/sbravyi/BivariateBicycleCodes):
# 8 layers per round. Neighbour index d of a check: X check -> 0,1,2 = left data via A1,A2,A3; 3,4,5 = right data via
# B1,B2,B3.  Z check -> 0,1,2 = left data via B1^T,B2^T,B3^T; 3,4,5 = right data via A1^T,A2^T,A3^T.
#   t=0      : RX on X checks            | CX data->Zcheck  (dir sZ[0])
#   t=1..5   : CX Xcheck->data (sX[t])    | CX data->Zcheck  (sZ[t])
#   t=6      : CX Xcheck->data (sX[6])    | M on Z checks
#   t=7      : MX on X checks, R on Z checks, all data idle
BRAVYI_SX = [None, 1, 4, 3, 5, 0, 2]
BRAVYI_SZ = [3, 5, 0, 1, 2, 4, None]


def bravyi_neighbours(name: str):
    """nbX[j][d], nbZ[i][d]: data-qubit index (left 0..lm-1, right lm..2lm-1) of neighbour d."""
    from .codes import bb_monomials
    A, B = bb_monomials(name)
    lm = A[0].shape[0]
    nbX = np.zeros((lm, 6), dtype=np.int64)
    nbZ = np.zeros((lm, 6), dtype=np.int64)
    for j in range(lm):
        for d in range(3):
            nbX[j, d] = int(np.nonzero(A[d][j, :])[0][0])
            nbX[j, 3 + d] = lm + int(np.nonzero(B[d][j, :])[0][0])
            nbZ[j, d] = int(np.nonzero(B[d][:, j])[0][0])            # B_d^T[j, :] = B_d[:, j]
            nbZ[j, 3 + d] = lm + int(np.nonzero(A[d][:, j])[0][0])   # A_d^T[j, :]
    return nbX, nbZ


def check_bravyi_schedule(code: CSSCode, nbX, nbZ, sX=BRAVYI_SX, sZ=BRAVYI_SZ):
    """Pure-numpy validity check of the interleaved schedule (raises AssertionError with a reason):
    (1) neighbour tables reproduce hx, hz; (2) layers t=1..5 touch every data qubit at most once;
    (3) for every (X check j, Z check i) the number of shared data qubits that X touches BEFORE Z is even
        (this makes every Z and X measurement equal to the stabilizer value; Stim re-checks it via the DEM)."""
    n = code.n
    lm = n // 2
    hx = np.zeros_like(code.hx)
    hz = np.zeros_like(code.hz)
    for j in range(lm):
        hx[j, nbX[j]] = 1
        hz[j, nbZ[j]] = 1
    assert np.array_equal(hx, code.hx), "X neighbour table does not reproduce hx"
    assert np.array_equal(hz, code.hz), "Z neighbour table does not reproduce hz"
    tX = -np.ones((lm, n), dtype=np.int64)
    tZ = -np.ones((lm, n), dtype=np.int64)
    for t in range(7):
        busy = []
        if sX[t] is not None:
            tX[np.arange(lm), nbX[:, sX[t]]] = t
            busy += nbX[:, sX[t]].tolist()
        if sZ[t] is not None:
            tZ[np.arange(lm), nbZ[:, sZ[t]]] = t
            busy += nbZ[:, sZ[t]].tolist()
        assert len(busy) == len(set(busy)), f"layer {t}: a data qubit is used twice"
    for j in range(lm):
        for i in range(lm):
            shared = (tX[j] >= 0) & (tZ[i] >= 0)
            nb = int(((tX[j] < tZ[i]) & shared).sum())
            assert nb % 2 == 0, f"X check {j} / Z check {i}: odd number ({nb}) of X-before-Z shared qubits"
    return True


def bb_memory_z_bravyi(name: str, code: CSSCode, p: float, rounds: int):
    """Memory-Z experiment with the depth-8 Bravyi cycle and the Bravyi noise model (uniform p):
    DEPOLARIZE2 after each CNOT, DEPOLARIZE1 on idle DATA qubits in every layer, flip after every reset,
    flip before every measurement. Data init |0> (+X flip), final Z readout of the data.
    Qubits: data 0..n-1, X checks n..n+mx-1, Z checks n+mx..n+mx+mz-1.
    Detectors: Z checks every round (round 0 vs |0>), X checks from round 1, final data vs last Z round."""
    import stim
    nbX, nbZ = bravyi_neighbours(name)
    check_bravyi_schedule(code, nbX, nbZ)
    n, mx, mz = code.n, code.hx.shape[0], code.hz.shape[0]
    data = list(range(n))
    xq = [n + i for i in range(mx)]
    zq = [n + mx + i for i in range(mz)]
    c = stim.Circuit()
    meas_index = {}
    counter = [0]

    def measure(qubits, keys, basis):
        c.append("X_ERROR" if basis == "z" else "Z_ERROR", qubits, p)
        c.append("M" if basis == "z" else "MX", qubits)
        for k in keys:
            meas_index[k] = counter[0]
            counter[0] += 1

    def rec(key):
        return stim.target_rec(meas_index[key] - counter[0])

    def cx(pairs):
        flat = [q for pr in pairs for q in pr]
        c.append("CX", flat)
        c.append("DEPOLARIZE2", flat, p)
        busy = set(flat)
        idle = [q for q in data if q not in busy]
        if idle:
            c.append("DEPOLARIZE1", idle, p)

    c.append("R", data + zq)
    c.append("X_ERROR", data + zq, p)
    c.append("TICK")
    sX, sZ = BRAVYI_SX, BRAVYI_SZ
    for r in range(rounds):
        # t = 0
        c.append("RX", xq)
        c.append("Z_ERROR", xq, p)
        cx([(int(nbZ[i, sZ[0]]), zq[i]) for i in range(mz)])
        c.append("TICK")
        # t = 1..5
        for t in range(1, 6):
            cx([(xq[j], int(nbX[j, sX[t]])) for j in range(mx)] + [(int(nbZ[i, sZ[t]]), zq[i]) for i in range(mz)])
            c.append("TICK")
        # t = 6
        measure(zq, [("z", i, r) for i in range(mz)], "z")
        cx([(xq[j], int(nbX[j, sX[6]])) for j in range(mx)])
        c.append("TICK")
        # t = 7
        measure(xq, [("x", j, r) for j in range(mx)], "x")
        if r < rounds - 1:
            c.append("R", zq)
            c.append("X_ERROR", zq, p)
        c.append("DEPOLARIZE1", data, p)
        c.append("TICK")
        for i in range(mz):
            if r == 0:
                c.append("DETECTOR", [rec(("z", i, r))], [i, r, 0])
            else:
                c.append("DETECTOR", [rec(("z", i, r)), rec(("z", i, r - 1))], [i, r, 0])
        if r > 0:
            for j in range(mx):
                c.append("DETECTOR", [rec(("x", j, r)), rec(("x", j, r - 1))], [j, r, 1])
    measure(data, [("d", q, 0) for q in data], "z")
    for i in range(mz):
        targets = [rec(("d", q, 0)) for q in np.nonzero(code.hz[i])[0]] + [rec(("z", i, rounds - 1))]
        c.append("DETECTOR", targets, [i, rounds, 0])
    for a in range(code.k):
        c.append("OBSERVABLE_INCLUDE", [rec(("d", q, 0)) for q in np.nonzero(code.lz[a])[0]], a)
    return c
