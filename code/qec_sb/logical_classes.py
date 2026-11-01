"""Pure binary-symplectic logical-class helpers, independent of any oracle."""

import numpy as np


def logical_class(code, ex, ez) -> int:
    ex = np.asarray(ex, dtype=np.uint8)
    ez = np.asarray(ez, dtype=np.uint8)
    lx = int((ez @ code.logical_x_support) & 1)
    lz = int((ex @ code.logical_z_support) & 1)
    return (lx << 1) | lz


def logical_class_seed(code, ex0, ez0, target_class):
    """Copy a feasible correction into the requested class using logical strings."""
    ex = np.asarray(ex0, dtype=np.uint8).copy()
    ez = np.asarray(ez0, dtype=np.uint8).copy()
    target_lx = (int(target_class) >> 1) & 1
    target_lz = int(target_class) & 1
    base_lx = int((ez @ code.logical_x_support) & 1)
    base_lz = int((ex @ code.logical_z_support) & 1)
    if base_lx != target_lx:
        ez ^= code.logical_z_support
    if base_lz != target_lz:
        ex ^= code.logical_x_support
    assert logical_class(code, ex, ez) == target_class
    return ex, ez
