from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True)
class RotatedSurfaceCode:
    distance: int
    Hx: np.ndarray  # X checks: detect Z component
    Hz: np.ndarray  # Z checks: detect X component
    logical_x_support: np.ndarray
    logical_z_support: np.ndarray
    data_coords: tuple

    @property
    def n(self) -> int:
        return self.distance * self.distance


def build_rotated_surface_code(distance: int) -> RotatedSurfaceCode:
    """Build the static rotated planar CSS geometry used by Stim's generator.

    This is a code-capacity parity-check model, not a circuit-level simulator.
    Coordinates/order mirror Stim's rotated surface-code generator.
    """
    d = int(distance)
    if d < 3 or d % 2 == 0:
        raise ValueError("Use an odd distance >= 3 for this prototype.")

    data = [(2 * ix + 1, 2 * iy + 1) for ix in range(d) for iy in range(d)]
    data_index = {q: i for i, q in enumerate(data)}
    x_meas, z_meas = [], []
    for x in range(d + 1):
        for y in range(d + 1):
            on_x_boundary = x == 0 or x == d
            on_y_boundary = y == 0 or y == d
            parity = (x % 2) != (y % 2)
            if on_x_boundary and parity:
                continue
            if on_y_boundary and not parity:
                continue
            (x_meas if parity else z_meas).append((2 * x, 2 * y))

    x_order = [(1, 1), (-1, 1), (1, -1), (-1, -1)]
    z_order = [(1, 1), (1, -1), (-1, 1), (-1, -1)]

    def support_matrix(measure_coords, order):
        rows = []
        for m in measure_coords:
            row = np.zeros(d * d, dtype=np.uint8)
            for dx, dy in order:
                q = (m[0] + dx, m[1] + dy)
                if q in data_index:
                    row[data_index[q]] = 1
            rows.append(row)
        return np.asarray(rows, dtype=np.uint8)

    Hx = support_matrix(x_meas, x_order)
    Hz = support_matrix(z_meas, z_order)
    logical_x = np.asarray([q[0] == 1 for q in data], dtype=np.uint8)
    logical_z = np.asarray([q[1] == 1 for q in data], dtype=np.uint8)
    return RotatedSurfaceCode(d, Hx, Hz, logical_x, logical_z, tuple(data))


def sample_depolarizing(code: RotatedSurfaceCode, p: float, shots: int, seed: int):
    """Sample code-capacity depolarizing Pauli errors.

    digit convention: 0=I, 1=X, 2=Z, 3=Y.
    Returns dictionaries containing X/Z components, syndromes and logical flips.
    """
    rng = np.random.default_rng(seed)
    probs = [1 - p, p / 3, p / 3, p / 3]
    pauli = rng.choice(4, size=(shots, code.n), p=probs)
    ex = ((pauli == 1) | (pauli == 3)).astype(np.uint8)
    ez = ((pauli == 2) | (pauli == 3)).astype(np.uint8)
    sx = (ex @ code.Hz.T) & 1
    sz = (ez @ code.Hx.T) & 1
    logical_x_flip = (ez @ code.logical_x_support) & 1
    logical_z_flip = (ex @ code.logical_z_support) & 1
    logical = np.stack([logical_x_flip, logical_z_flip], axis=1).astype(np.uint8)
    return {"pauli": pauli, "ex": ex, "ez": ez, "sx": sx, "sz": sz, "logical": logical}
