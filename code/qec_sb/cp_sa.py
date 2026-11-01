"""Constraint-preserving Metropolis annealing in stabilizer-generator bits."""

from dataclasses import dataclass
from math import exp, isfinite
from numbers import Integral, Real
from time import perf_counter, process_time

import numpy as np

from . import cp_sb as _cp_sb
from .cp_sb import StabilizerHamiltonian


@dataclass(frozen=True)
class CPSAParams:
    """Annealing controls; ``chains`` is the number for each restart mode."""

    chains: int = 8
    sweeps: int = 100
    temperature_start: float = 2.0
    temperature_end: float = 0.05
    restart_modes: tuple[str, ...] = ("warm", "random")

    def __post_init__(self):
        for name in ("chains", "sweeps"):
            value = getattr(self, name)
            if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral) or value <= 0:
                raise ValueError(f"{name} must be a positive integer")
        for name in ("temperature_start", "temperature_end"):
            value = getattr(self, name)
            if (isinstance(value, (bool, np.bool_)) or not isinstance(value, Real)
                    or not isfinite(value) or value <= 0):
                raise ValueError(f"{name} must be finite and strictly positive")
        if self.temperature_end > self.temperature_start:
            raise ValueError("temperature_end must not exceed temperature_start")
        if (not isinstance(self.restart_modes, tuple) or not self.restart_modes
                or any(mode not in ("warm", "random") for mode in self.restart_modes)):
            raise ValueError("restart_modes must be a nonempty tuple containing 'warm' or 'random'")
        if len(set(self.restart_modes)) != len(self.restart_modes):
            raise ValueError("restart_modes must not contain duplicates")


@dataclass(frozen=True)
class CPSAResult:
    """Annealing outputs with shallow freezing of field bindings only.

    Array contents and the dictionaries in ``trace`` remain mutable. Callers
    needing independent editable data must copy the arrays and trace dictionaries
    before modifying their contents; the frozen dataclass does not protect them.
    """

    terminal_bits: np.ndarray
    terminal_weights: np.ndarray
    best_ever_bits: np.ndarray
    best_ever_weights: np.ndarray
    proposals: int
    accepted: int
    energy_evaluations: int
    wall_time_s: float
    process_time_s: float
    trace: tuple[dict, ...]


def _single_bits(ham, bits, name):
    bits = np.asarray(bits)
    if bits.shape != (ham.n_spins,):
        raise ValueError(f"{name} must have shape [n_spins]")
    if not np.all((bits == 0) | (bits == 1)):
        raise ValueError(f"{name} must contain only binary values 0 or 1")
    return bits.astype(np.uint8, copy=True)


def _initialize_bits(ham, params, rng, warm_bits=None):
    """Return mode-major states; omitted warm bits conservatively mean u=0."""
    warm = (np.zeros(ham.n_spins, dtype=np.uint8) if warm_bits is None
            else _single_bits(ham, warm_bits, "warm_bits"))
    states = []
    for mode in params.restart_modes:
        if mode == "warm" or ham.n_spins == 0:
            states.append(np.tile(warm, (params.chains, 1)))
        else:
            states.append(rng.integers(0, 2, size=(params.chains, ham.n_spins), dtype=np.uint8))
    return np.concatenate(states, axis=0)


def _temperature_schedule(params):
    """Geometric cooling with both endpoints; a single sweep uses the start."""
    start, end = float(params.temperature_start), float(params.temperature_end)
    if params.sweeps == 1:
        return np.asarray([start])
    if start == end:
        return np.full(params.sweeps, start)
    log_temperatures = np.linspace(np.log(start), np.log(end), params.sweeps)
    # Clamp roundoff at the representable limits, then restore exact endpoints.
    with np.errstate(over="ignore", under="ignore"):
        temperatures = np.exp(log_temperatures)
    np.clip(temperatures, end, start, out=temperatures)
    temperatures[0], temperatures[-1] = start, end
    return temperatures


def _physical_weight(physical, n_data):
    return int(np.count_nonzero(physical[:n_data] | physical[n_data:]))


def proposal_delta_bits(ham: StabilizerHamiltonian, u, spin_index):
    """Exact E(u XOR e_j)-E(u) from XORing physical generator j, without mutation."""
    bits = _single_bits(ham, u, "u")
    if (isinstance(spin_index, (bool, np.bool_)) or not isinstance(spin_index, Integral)
            or not 0 <= spin_index < ham.n_spins):
        raise ValueError("spin_index must be an integer in [0, n_spins)")
    physical = ham.physical_from_bits(bits)
    candidate = physical ^ ham.generators[spin_index]
    return _physical_weight(candidate, ham.n_data) - _physical_weight(physical, ham.n_data)


def run_constraint_preserving_sa(
        ham: StabilizerHamiltonian, params: CPSAParams, seed: int,
        warm_bits=None, return_trace=False):
    """Anneal within one syndrome/logical class using exact physical Pauli weight.

    Each restart mode contributes ``params.chains`` rows, in the requested mode
    order. Warm rows copy ``warm_bits`` or use u=0 (the Hamiltonian seed) if absent.
    Random rows are independent seeded binary states. Every sweep visits a fresh
    seeded random permutation of all spins for each chain. Temperatures cool
    geometrically, including both positive endpoints when there are >=2 sweeps;
    one sweep uses ``temperature_start``.

    Non-increasing moves are always accepted. Increasing moves are accepted with
    probability exp(-delta/T). Best-ever states include initialization and retain
    the first state at a tied best weight; terminal states are recorded separately.

    ``proposals`` counts visited spin flips, ``accepted`` counts committed flips,
    and ``energy_evaluations`` counts one initial weight per chain plus one exact
    candidate weight per proposal. Current weights are cached. Trace summaries
    use cached weights, draw no randomness, and add no energy evaluations. With
    zero spins, chains remain at the seed with no proposals or random draws.
    Timing includes initialization, annealing, and optional trace collection;
    all other result fields are reproducible for the same inputs and seed.
    """
    wall_start, process_start = perf_counter(), process_time()
    rng = np.random.default_rng(seed)
    bits = _initialize_bits(ham, params, rng, warm_bits)
    physical = ham.physical_from_bits(bits)
    weights = np.count_nonzero(
        physical[:, :ham.n_data] | physical[:, ham.n_data:], axis=1
    ).astype(int)
    best_bits, best_weights = bits.copy(), weights.copy()
    proposals = accepted = 0
    energy_evaluations = len(bits)
    trace = []

    # dE of a proposed flip: "sparse" (default) reads only the <= K qubits of the generator; "dense" (original)
    # XORs and recounts the full 2n-bit error. Same integers, same random draws -> identical results.
    sparse = _cp_sb.DELTA_MODE != "dense"
    sup_lists = _cp_sb.generator_support(ham).lists if (sparse and ham.n_spins) else None
    n_data = ham.n_data

    for sweep, temperature in enumerate(_temperature_schedule(params), start=1):
        if ham.n_spins and sparse:
            T = float(temperature)
            for chain in range(len(bits)):
                exl = physical[chain, :n_data].tolist()
                ezl = physical[chain, n_data:].tolist()
                w = int(weights[chain])
                brow = bits[chain]
                for spin in rng.permutation(ham.n_spins):
                    sup = sup_lists[spin]
                    delta = 0
                    for q, fx, fz in sup:
                        a, b = exl[q], ezl[q]
                        delta += ((a ^ fx) | (b ^ fz)) - (a | b)
                    proposals += 1
                    energy_evaluations += 1
                    if delta <= 0 or rng.random() < exp(-delta / T):
                        for q, fx, fz in sup:
                            exl[q] ^= fx
                            ezl[q] ^= fz
                        brow[spin] ^= 1
                        w += delta
                        accepted += 1
                        if w < best_weights[chain]:
                            best_bits[chain] = brow
                            best_weights[chain] = w
                physical[chain, :n_data] = exl
                physical[chain, n_data:] = ezl
                weights[chain] = w
        elif ham.n_spins:
            for chain in range(len(bits)):
                for spin in rng.permutation(ham.n_spins):
                    candidate = physical[chain] ^ ham.generators[spin]
                    candidate_weight = _physical_weight(candidate, ham.n_data)
                    delta = candidate_weight - int(weights[chain])
                    proposals += 1
                    energy_evaluations += 1
                    if delta <= 0 or rng.random() < exp(-delta / float(temperature)):
                        bits[chain, spin] ^= 1
                        physical[chain] = candidate
                        weights[chain] = candidate_weight
                        accepted += 1
                        if candidate_weight < best_weights[chain]:
                            best_bits[chain] = bits[chain]
                            best_weights[chain] = candidate_weight
        if return_trace:
            trace.append({
                "sweep": sweep,
                "temperature": float(temperature),
                "best_weight": float(weights.min()),
                "mean_weight": float(weights.mean()),
                "best_ever_weight": float(best_weights.min()),
                "feasible_fraction": 1.0,
            })

    return CPSAResult(
        terminal_bits=bits,
        terminal_weights=weights,
        best_ever_bits=best_bits,
        best_ever_weights=best_weights,
        proposals=proposals,
        accepted=accepted,
        energy_evaluations=energy_evaluations,
        wall_time_s=perf_counter() - wall_start,
        process_time_s=process_time() - process_start,
        trace=tuple(trace),
    )
