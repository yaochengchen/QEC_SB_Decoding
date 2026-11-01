from dataclasses import dataclass
import numpy as np

@dataclass
class SBParams:
    replicas: int = 64
    steps: int = 800
    dt: float = 0.5
    warm_strength: float = 0.08
    gb_control_A: float = 0.2


def _coupling_2021(J):
    N = J.shape[0]
    mask = ~np.eye(N, dtype=bool)
    jrms = np.sqrt(np.sum(J[mask] ** 2) / max(N * (N - 1), 1))
    return 0.5 / max(jrms * np.sqrt(N), 1e-12)


def _gb_coupling_and_dt(J, requested_dt):
    vals = np.linalg.eigvalsh(J)
    lmax = max(float(vals[-1]), 1e-12)
    lmin = float(vals[0])
    c = 1.0 / lmax
    if requested_dt is None:
        # Goto-Hidaka-Tatsumura 2026 Appendix-A spectral stability heuristic.
        requested_dt = 1.25 * np.sqrt(2.0 / max(1.0 - lmin / lmax, 1e-12))
    return c, float(requested_dt)


def run_sb(J, variant: str, params: SBParams, seed: int, warm_spin=None,
           return_trace=False, trace_metrics=None):
    """Run bSB, dSB, or generalized ballistic SB (GbSB).

    bSB/dSB implement the 2021 symplectic-Euler hard-wall equations.
    GbSB implements the 2026 individualized bifurcation-parameter update:
      p_i(m+1)=p_i(m)-[1-A x_i(m)^2] p_i(m)/(M-m).
    When ``return_trace`` is true, return ``(spins, trace)`` with the best and
    mean Ising objective after every iteration. An optional read-only
    ``trace_metrics(spins)`` callback may add domain-specific metrics such as
    syndrome-feasible fraction. Trace collection does not alter random draws or
    the dynamics.
    """
    variant = variant.lower()
    if variant not in {"bsb", "dsb", "gbsb"}:
        raise ValueError("variant must be bsb, dsb, or gbsb")
    rng = np.random.default_rng(seed)
    R, N = params.replicas, J.shape[0]
    x = rng.uniform(-0.1, 0.1, size=(R, N))
    y = np.zeros((R, N), dtype=float)
    if warm_spin is not None:
        ws = np.asarray(warm_spin, dtype=float)
        x = params.warm_strength * ws[None, :] + rng.normal(0, 0.01, size=(R, N))

    if variant == "gbsb":
        c, dt = _gb_coupling_and_dt(J, params.dt)
        p_ind = np.ones((R, N), dtype=float)
    else:
        c, dt = _coupling_2021(J), float(params.dt)

    M = params.steps
    trace = [] if return_trace else None
    for m in range(M):
        if variant == "gbsb":
            denom = max(M - m, 1)
            p_ind -= (1.0 - params.gb_control_A * x * x) * p_ind / denom
            interaction = x @ J.T
            y += (-p_ind * x + c * interaction) * dt
        else:
            a = (m + 1) / M  # a(t): 0 -> a0=1
            restoring = -(1.0 - a) * x
            source = np.sign(x) if variant == "dsb" else x
            interaction = source @ J.T
            y += (restoring + c * interaction) * dt
        x += y * dt
        hit = np.abs(x) > 1.0
        x = np.where(hit, np.sign(x), x)
        y = np.where(hit, 0.0, y)

        if return_trace:
            trace_spins = np.sign(x)
            trace_spins[trace_spins == 0] = 1.0
            energies = -0.5 * np.einsum(
                "ri,ij,rj->r", trace_spins, J, trace_spins
            )
            row = {
                "iteration": m + 1,
                "best_energy": float(np.min(energies)),
                "mean_energy": float(np.mean(energies)),
            }
            if trace_metrics is not None:
                extra = dict(trace_metrics(trace_spins.astype(np.int8)))
                overlap = set(row).intersection(extra)
                if overlap:
                    raise ValueError(f"trace metric overrides reserved fields: {sorted(overlap)}")
                row.update(extra)
            trace.append(row)

    spin = np.sign(x)
    spin[spin == 0] = 1
    spins = spin.astype(np.int8)
    return (spins, trace) if return_trace else spins


def binary_from_augmented_spins(spins):
    """Gauge-fix the auxiliary spin and return original binary variables."""
    effective = spins[:, :-1] * spins[:, -1:]
    return ((1 - effective) // 2).astype(np.uint8)
