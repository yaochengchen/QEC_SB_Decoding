from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from numbers import Integral, Real
from time import perf_counter, process_time
from typing import Dict, Mapping, Tuple
import numpy as np
import networkx as nx

from .rotated_code import RotatedSurfaceCode
from .logical_classes import logical_class as _logical_class, logical_class_seed as _class_seed
from .cp_sb import build_stabilizer_hamiltonian, run_constraint_preserving_sb, CPSBParams


class SolverOutputValidationError(ValueError):
    """A malformed returned solver value, distinct from execution exceptions."""

    def __init__(self, code, message):
        self.code = code
        super().__init__(message)


def validate_solver_number(value, label):
    """Validate a returned real, finite, nonnegative scalar before arithmetic."""
    valid = not isinstance(value, (bool, np.bool_)) and isinstance(value, Real)
    try:
        valid = valid and isfinite(value) and value >= 0
    except OverflowError:
        valid = False
    if not valid:
        raise SolverOutputValidationError("invalid_solver_number", f"{label} must be real, finite and nonnegative")
    return float(value)


def validate_solver_counter(value, label, *, expected=None):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral) or value < 0:
        raise SolverOutputValidationError("invalid_solver_counter", f"{label} must be a nonnegative integer, not bool")
    if expected is not None and value != expected:
        raise SolverOutputValidationError("inconsistent_solver_counter", f"{label} disagrees with the requested native work")
    return int(value)


def validate_solver_array(value, shape, label, *, allowed=None):
    try:
        array = np.asarray(value)
    except (TypeError, ValueError) as exc:
        raise SolverOutputValidationError("invalid_solver_array", f"{label} must be a real numeric array") from exc
    if array.dtype.kind not in "iuf":
        raise SolverOutputValidationError("invalid_solver_array", f"{label} must be a real numeric array")
    if array.shape != shape:
        raise SolverOutputValidationError("invalid_solver_array", f"{label} must have shape {shape}")
    if not np.all(np.isfinite(array)):
        raise SolverOutputValidationError("invalid_solver_array", f"{label} must all be finite")
    if allowed is not None and not np.all(np.isin(array, allowed)):
        raise SolverOutputValidationError("invalid_solver_array", f"{label} contains values outside {allowed}")
    return array


def validate_solver_mapping(value, required, label, *, exact=False):
    if not isinstance(value, Mapping):
        raise SolverOutputValidationError("invalid_solver_mapping", f"{label} must be a mapping")
    missing = set(required) - set(value)
    if missing:
        raise SolverOutputValidationError("invalid_solver_mapping", f"{label} is missing {', '.join(map(str, sorted(missing)))}")
    if exact and set(value) != set(required):
        raise SolverOutputValidationError("invalid_solver_mapping", f"{label} has unexpected keys")
    return value


def _validate_sb_result(output, ham, params, diagnostics):
    if diagnostics:
        if not isinstance(output, (tuple, list)) or len(output) != 2:
            raise SolverOutputValidationError("invalid_sb_result", "SB result must be a (spins, trace) pair")
        spins, trace = output
    else:
        spins, trace = output, []
    spins = validate_solver_array(spins, (params.replicas, ham.n_spins), "SB solver spins", allowed=(-1, 1))
    if diagnostics:
        if not isinstance(trace, (tuple, list)) or (len(trace) != params.steps and not (ham.n_spins == 0 and len(trace) == 0)):
            raise SolverOutputValidationError("invalid_sb_trace", "SB trace must contain one mapping per requested iteration")
        for index, row in enumerate(trace, start=1):
            validate_solver_mapping(row, ("iteration", "best_weight", "mean_weight", "feasible_fraction"), "SB trace row")
            validate_solver_counter(row["iteration"], "SB trace iteration", expected=index)
            best = validate_solver_number(row["best_weight"], "SB trace best_weight")
            mean = validate_solver_number(row["mean_weight"], "SB trace mean_weight")
            feasible = validate_solver_number(row["feasible_fraction"], "SB trace feasible_fraction")
            if best > mean or mean > ham.n_data or feasible > 1:
                raise SolverOutputValidationError("invalid_sb_trace", "SB trace weight/feasible_fraction is outside physical bounds")
    return spins, trace


def gf2_solve(A: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Return one solution x of A x=b over GF(2), with free vars set to zero."""
    A = np.asarray(A, dtype=np.uint8).copy()
    b = np.asarray(b, dtype=np.uint8).copy().reshape(-1)
    m, n = A.shape
    aug = np.concatenate([A, b[:, None]], axis=1)
    pivots = []
    r = 0
    for c in range(n):
        rows = np.flatnonzero(aug[r:, c])
        if len(rows) == 0:
            continue
        p = r + int(rows[0])
        if p != r:
            aug[[r, p]] = aug[[p, r]]
        for rr in range(m):
            if rr != r and aug[rr, c]:
                aug[rr] ^= aug[r]
        pivots.append(c)
        r += 1
        if r == m:
            break
    for rr in range(r, m):
        if not aug[rr, :n].any() and aug[rr, n]:
            raise ValueError("Inconsistent GF(2) system")
    x = np.zeros(n, dtype=np.uint8)
    for rr, c in enumerate(pivots):
        x[c] = aug[rr, n]
    if not np.array_equal((A @ x) & 1, b):
        raise RuntimeError("GF(2) solver verification failed")
    return x


def feasible_correction(code: RotatedSurfaceCode, sx, sz):
    ex = gf2_solve(code.Hz, np.asarray(sx, dtype=np.uint8))
    ez = gf2_solve(code.Hx, np.asarray(sz, dtype=np.uint8))
    return ex, ez


class _ComponentMWPM:
    """Code-capacity unit-weight matching decoder for one CSS component.

    H has check rows and data-qubit columns, each column of weight 1 or 2.
    Weight-1 columns are boundaries. This implementation is deliberately small,
    transparent, and depends only on NetworkX; it is used as the common MWPM
    baseline and as a high-quality feasible seed when PyMatching is unavailable.
    """
    def __init__(self, H: np.ndarray):
        H = np.asarray(H, dtype=np.uint8)
        self.H = H
        self.n_checks, self.n_qubits = H.shape
        self.graph = nx.Graph()
        self.graph.add_nodes_from(range(self.n_checks))
        self.boundary_edges: Dict[int, list[int]] = {i: [] for i in range(self.n_checks)}
        self.edge_qubit: Dict[Tuple[int,int], int] = {}
        for q in range(self.n_qubits):
            checks = np.flatnonzero(H[:, q]).tolist()
            if len(checks) == 2:
                a, b = int(checks[0]), int(checks[1])
                key = (min(a,b), max(a,b))
                # Rotated surface-code geometry has no parallel check edges.
                self.graph.add_edge(a, b, weight=1)
                self.edge_qubit[key] = q
            elif len(checks) == 1:
                self.boundary_edges[int(checks[0])].append(q)
            elif len(checks) != 0:
                raise ValueError("Parity matrix is not graph-like (column weight >2)")
        # All-pairs paths on the small check graph.
        self.paths = dict(nx.all_pairs_shortest_path(self.graph))
        self.dists = dict(nx.all_pairs_shortest_path_length(self.graph))
        boundary_checks = [c for c, qs in self.boundary_edges.items() if qs]
        self.best_boundary = {}
        for s in range(self.n_checks):
            best = None
            for b in boundary_checks:
                if b not in self.dists.get(s, {}):
                    continue
                cost = self.dists[s][b] + 1
                if best is None or cost < best[0] or (cost == best[0] and b < best[1]):
                    best = (cost, b, self.boundary_edges[b][0])
            if best is None:
                raise ValueError("No boundary reachable from check graph")
            self.best_boundary[s] = best

    def _path_qubits(self, a: int, b: int):
        nodes = self.paths[a][b]
        qs = []
        for u, v in zip(nodes[:-1], nodes[1:]):
            qs.append(self.edge_qubit[(min(u,v), max(u,v))])
        return qs

    def decode(self, syndrome: np.ndarray) -> np.ndarray:
        syndrome = np.asarray(syndrome, dtype=np.uint8)
        defects = np.flatnonzero(syndrome).tolist()
        corr = np.zeros(self.n_qubits, dtype=np.uint8)
        k = len(defects)
        if k == 0:
            return corr
        G = nx.Graph()
        D = [("d", i) for i in range(k)]
        B = [("b", i) for i in range(k)]
        G.add_nodes_from(D + B)
        for i in range(k):
            a = defects[i]
            for j in range(i + 1, k):
                b = defects[j]
                if b not in self.dists.get(a, {}):
                    continue
                G.add_edge(D[i], D[j], weight=int(self.dists[a][b]))
            G.add_edge(D[i], B[i], weight=int(self.best_boundary[a][0]))
        for i in range(k):
            for j in range(i + 1, k):
                G.add_edge(B[i], B[j], weight=0)
        matching = nx.min_weight_matching(G, weight="weight")
        if len(matching) * 2 != len(G):
            raise RuntimeError("MWPM did not produce a perfect matching")
        for u, v in matching:
            if u[0] == "b" and v[0] == "b":
                continue
            if u[0] == "b":
                u, v = v, u
            if u[0] != "d":
                continue
            i = int(u[1]); a = defects[i]
            if v[0] == "d":
                b = defects[int(v[1])]
                for q in self._path_qubits(a, b):
                    corr[q] ^= 1
            else:
                _, bc, bq = self.best_boundary[a]
                for q in self._path_qubits(a, bc):
                    corr[q] ^= 1
                corr[bq] ^= 1
        if not np.array_equal((self.H @ corr) & 1, syndrome):
            raise RuntimeError("MWPM correction has wrong syndrome")
        return corr


class CSSMWPMDecoder:
    def __init__(self, code: RotatedSurfaceCode):
        # X error component is checked by Hz; Z error component by Hx.
        self.code = code
        self.xdec = _ComponentMWPM(code.Hz)
        self.zdec = _ComponentMWPM(code.Hx)

    def decode(self, sx, sz):
        ex = self.xdec.decode(np.asarray(sx, dtype=np.uint8))
        ez = self.zdec.decode(np.asarray(sz, dtype=np.uint8))
        return ex, ez


@dataclass
class ScalableCPConfig:
    replicas: int = 12
    steps: int = 160
    dt: float = 0.45
    coupling_scale: float = 0.70
    warm_strength: float = 0.08
    random_restarts: int = 1
    warm_restarts: int = 1
    coupling_multipliers: tuple = (0.85, 1.0, 1.15)


def collect_cp_class_pools_no_oracle(code, ex0, ez0, base_seed: int,
                                     cfg: ScalableCPConfig, *, hams=None,
                                     seed_for_run=None, diagnostics=False):
    """Collect terminal SB states plus u=0 with the established restart ordering.

    Optional prebuilt ``hams`` let a comparison share the exact model instances.
    ``seed_for_run(class_id, run_id)`` overrides only the legacy seed arithmetic.
    Diagnostics retain raw terminal trajectories and post-step trace minima;
    they never add intermediate states to the prediction pool.
    ``force_evaluations`` counts timestep batch queries (table lookups at d=3),
    ``force_vector_evaluations`` counts replica queries, and table construction
    is counted separately. None of these units equates to SA spin proposals.
    """
    pools = {}
    build_hams = hams is None
    hams = {} if build_hams else hams
    stats = {}
    for c in range(4):
        if build_hams:
            cex, cez = _class_seed(code, ex0, ez0, c)
            hams[c] = build_stabilizer_hamiltonian(code, cex, cez)
        ham = hams[c]
        raw = [np.zeros((1, ham.n_spins), dtype=np.uint8)]
        run_id = 0
        wall_s = process_s = 0.0
        restart_metrics = []
        best_observed = None
        for mult in cfg.coupling_multipliers:
            par = CPSBParams(
                replicas=cfg.replicas,
                steps=cfg.steps,
                dt=cfg.dt,
                coupling_scale=cfg.coupling_scale * float(mult),
                warm_strength=cfg.warm_strength,
            )
            for warm, repeats in ((True, cfg.warm_restarts), (False, cfg.random_restarts)):
                for _ in range(int(repeats)):
                    run_seed = (int(base_seed + c*100003 + run_id*997) if seed_for_run is None
                                else seed_for_run(c, run_id))
                    wall_start, process_start = perf_counter(), process_time()
                    out = run_constraint_preserving_sb(
                        ham, par, run_seed,
                        warm_bits=np.zeros(ham.n_spins, dtype=np.uint8) if warm else None,
                        return_trace=diagnostics,
                    )
                    # Validate each measured interval before summing: a later
                    # positive interval must not cancel a negative one.
                    wall_elapsed = validate_solver_number(perf_counter() - wall_start, "SB restart wall_time_s")
                    process_elapsed = validate_solver_number(process_time() - process_start, "SB restart process_time_s")
                    native = {
                        "wall_time_s": wall_elapsed, "process_time_s": process_elapsed,
                        "force_evaluations": validate_solver_counter(cfg.steps, "SB restart force_evaluations"),
                        "force_vector_evaluations": validate_solver_counter(cfg.steps * cfg.replicas, "SB restart force_vector_evaluations"),
                        "force_table_precompute_vectors": (1 << ham.n_spins) if ham.n_spins <= 16 else 0,
                    }
                    s, trace = _validate_sb_result(out, ham, par, diagnostics)
                    restart_metrics.append(native)
                    wall_s += wall_elapsed
                    process_s += process_elapsed
                    if diagnostics:
                        observed = min((row["best_weight"] for row in trace), default=None)
                        if observed is not None:
                            best_observed = (observed if best_observed is None
                                             else min(best_observed, observed))
                    raw.append((s < 0).astype(np.uint8))
                    run_id += 1
        u = np.unique(np.vstack(raw), axis=0)
        phys = ham.physical_from_bits(u)
        wt = np.count_nonzero(phys[:, :code.n] | phys[:, code.n:], axis=1)
        order = np.argsort(wt, kind="stable")
        pools[c] = u[order]
        stats[c] = {
            "unique_cp_candidates": int(len(u)),
            "best_weight": int(wt[order[0]]),
            "mean_weight": float(wt.mean()),
        }
        if diagnostics:
            terminal = np.vstack(raw[1:]) if len(raw) > 1 else np.empty((0, ham.n_spins), np.uint8)
            unique_terminal = len(np.unique(terminal, axis=0))
            stats[c].update({
                "terminal_bits": terminal,
                "terminal_trajectories": len(terminal),
                "unique_terminal_candidates": unique_terminal,
                "terminal_duplicate_rate": 1 - unique_terminal / len(terminal) if len(terminal) else 0.0,
                "best_observed_step_weight": best_observed,
                "best_ever_bits": None,
                "best_ever_weights": None,
                "best_ever_scope": "SB trace minimum after updates; initialization is not observed",
                "proposals": None,
                "accepted": None,
                "energy_evaluations": None,
                "force_evaluations": run_id * cfg.steps,
                "force_vector_evaluations": run_id * cfg.steps * cfg.replicas,
                "force_table_precompute_vectors": run_id * (1 << ham.n_spins) if ham.n_spins <= 16 else 0,
                "wall_time_s": wall_s,
                "process_time_s": process_s,
                "restart_metrics": restart_metrics,
            })
    return pools, hams, stats


def cp_map_prediction(pools, hams):
    best = None
    for c in range(4):
        u = pools[c]
        phys = hams[c].physical_from_bits(u)
        n = hams[c].n_data
        wt = np.count_nonzero(phys[:, :n] | phys[:, n:], axis=1)
        w = int(wt.min())
        if best is None or w < best[0] or (w == best[0] and c < best[1]):
            best = (w, c)
    return int(best[1]), int(best[0])


def correction_class(code, ex, ez):
    return _logical_class(code, np.asarray(ex,dtype=np.uint8), np.asarray(ez,dtype=np.uint8))
