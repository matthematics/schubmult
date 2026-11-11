"""Experimental combinatorial models.

These implementations are exploratory and may be incomplete or incorrect.
Their classes are available here or from their individual submodules.
"""

_exports = {
    "ChuteMoveElement": "chute_move_element",
    "HeckePlactic": "hecke_plactic",
    "HPD": "hpd",
    "HPDTile": "hpd",
    "InversionsTableau": "inversions_tableau",
    "QuasiCrystalGraph": "quasi_crystal_graph",
    "QuasiCrystalGraphTensor": "quasi_crystal_graph",
    "RootTableau": "root_tableau",
}

__all__ = list(_exports)


def __getattr__(name):
    if name not in _exports:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    import importlib

    value = getattr(importlib.import_module(f"{__name__}.{_exports[name]}"), name)
    globals()[name] = value
    return value


def __dir__():
    return sorted(set(globals()) | set(__all__))
