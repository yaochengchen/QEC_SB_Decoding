from dataclasses import dataclass
import numpy as np
from .rotated_code import RotatedSurfaceCode


def syndrome_key(sx: np.ndarray, sz: np.ndarray) -> int:
    out = 0
    for bit in np.concatenate([sx, sz]):
        out = (out << 1) | int(bit)
    return out


def _logical_tuple(code, ex, ez):
    return (int((ez @ code.logical_x_support) & 1), int((ex @ code.logical_z_support) & 1))

@dataclass
class D3Oracle:
    min_pauli_weight: dict
    map_config: dict
    independent_x: dict
    independent_z: dict
    # key -> counts[logical_class(0..3), pauli_weight]
    class_weight_counts: dict

    def exact_mld_class(self, key: int, p: float):
        counts = self.class_weight_counts[key]
        r = (p / 3.0) / (1.0 - p)
        powers = r ** np.arange(counts.shape[1])
        masses = counts @ powers
        cls = int(np.argmax(masses))
        return np.asarray([cls >> 1, cls & 1], dtype=np.uint8), masses


def build_d3_oracle(code: RotatedSurfaceCode) -> D3Oracle:
    """Exact enumeration oracle for d=3 (4^9 = 262,144 Pauli errors)."""
    if code.distance != 3:
        raise ValueError("Exact enumeration oracle is intentionally limited to d=3.")
    n = code.n

    # Independent X/Z minimum-weight decoders (the uniform code-capacity MWPM objective).
    bit_ids = np.arange(1 << n, dtype=np.uint16)[:, None]
    bits = ((bit_ids >> np.arange(n)) & 1).astype(np.uint8)
    sx_all = (bits @ code.Hz.T) & 1
    sz_all = (bits @ code.Hx.T) & 1
    independent_x, independent_z = {}, {}
    for syn_all, store in ((sx_all, independent_x), (sz_all, independent_z)):
        for i, e in enumerate(bits):
            key = tuple(int(v) for v in syn_all[i])
            wt = int(e.sum())
            if key not in store or wt < store[key][0]:
                store[key] = (wt, e.copy())

    min_w, map_config, counts_by_key = {}, {}, {}
    total = 4 ** n
    powers4 = 4 ** np.arange(n, dtype=np.int64)
    chunk = 32768
    for start in range(0, total, chunk):
        ids = np.arange(start, min(total, start + chunk), dtype=np.int64)
        pauli = ((ids[:, None] // powers4) % 4).astype(np.uint8)
        ex = ((pauli == 1) | (pauli == 3)).astype(np.uint8)
        ez = ((pauli == 2) | (pauli == 3)).astype(np.uint8)
        sx = (ex @ code.Hz.T) & 1
        sz = (ez @ code.Hx.T) & 1
        wt = (pauli != 0).sum(axis=1)
        lx = (ez @ code.logical_x_support) & 1
        lz = (ex @ code.logical_z_support) & 1
        for i in range(len(ids)):
            key = syndrome_key(sx[i], sz[i])
            wi = int(wt[i])
            if key not in min_w or wi < min_w[key]:
                min_w[key] = wi
                map_config[key] = (ex[i].copy(), ez[i].copy())
            if key not in counts_by_key:
                counts_by_key[key] = np.zeros((4, n + 1), dtype=np.int64)
            cls = (int(lx[i]) << 1) | int(lz[i])
            counts_by_key[key][cls, wi] += 1

    return D3Oracle(min_w, map_config, independent_x, independent_z, counts_by_key)


def independent_decode(code: RotatedSurfaceCode, oracle: D3Oracle, sx, sz):
    ex = oracle.independent_x[tuple(int(v) for v in sx)][1].copy()
    ez = oracle.independent_z[tuple(int(v) for v in sz)][1].copy()
    logical = np.asarray(_logical_tuple(code, ex, ez), dtype=np.uint8)
    return ex, ez, logical
