"""QEC-aware simulated-bifurcation research prototype."""
__version__ = "0.2.0"

# Preserve the P4.1 package exports without loading its exact-oracle experiment
# dependencies when importing independent solver/predictor modules.
__all__ = ["ImportanceMassConfig", "P41Config", "adaptive_p41_prediction"]


def __getattr__(name):
    if name in __all__:
        from . import importance_mass
        value = getattr(importance_mass, name)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(set(globals()) | set(__all__))
