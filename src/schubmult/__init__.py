import importlib
import inspect
import pkgutil
import warnings
from typing import List

__version__ = "4.0.0dev"

"""
Top-level schubmult package initializer.

Dynamically import and re-export public classes defined in submodules
under this package. Imports are best-effort: modules that fail to import
are skipped so that a top-level `import schubmult` does not raise due to
optional heavy dependencies (e.g. pulp/mosek).
"""

__all__: List[str] = []
# map exported name -> module that originally defined it (for lazy re-importing)
_module_map = {}

# Walk submodules and import them where possible. Skip scripts and skip on errors.
for finder, modname, ispkg in pkgutil.walk_packages(__path__, prefix=__name__ + "."):
    if modname.startswith(__name__ + ".scripts"):
        continue
    try:
        mod = __import__(modname, fromlist=["*"])
    except Exception as e:
        # Don't fail import of the whole package if a submodule has heavy deps.
        warnings.warn(f"skipping import of {modname!r} due to error: {e!r}", RuntimeWarning)
        continue

    for attr_name, attr_val in vars(mod).items():
        if attr_name.startswith("_"):
            continue
        if inspect.isclass(attr_val) and getattr(attr_val, "__module__", None) == mod.__name__:
            globals()[attr_name] = attr_val
            __all__.append(attr_name)
            _module_map[attr_name] = mod.__name__

# Commonly exported helpers that we prefer to lazy-load
# we do not import them eagerly to avoid triggering heavy dependencies.
_lazy_exports = {
    "Sx": "schubmult.rings.schubert_ring",
    "uncode": "schubmult.schub_lib.perm_lib",
}
for name in _lazy_exports:
    if name not in __all__:
        __all__.append(name)

# finalize __all__
__all__ = tuple(sorted(set(__all__)))

def __getattr__(name: str):
    """
    Lazy import exported names that were not successfully imported at package
    initialization time. This allows `import schubmult` to succeed even if
    some optional dependencies (used by particular submodules) are missing.
    """
    if name in globals():
        return globals()[name]
    # known lazy exports
    if name in _lazy_exports:
        try:
            mod = importlib.import_module(_lazy_exports[name])
            val = getattr(mod, name)
            globals()[name] = val
            return val
        except Exception as e:
            raise AttributeError(f"cannot import {name!r} from {_lazy_exports[name]!r}: {e}") from e
    # names discovered during the initial scan
    if name in _module_map:
        try:
            mod = importlib.import_module(_module_map[name])
            val = getattr(mod, name)
            globals()[name] = val
            return val
        except Exception as e:
            raise AttributeError(f"cannot import {name!r} from {_module_map[name]!r}: {e}") from e
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

def __dir__():
    # include lazily-exported names in dir()
    return sorted(set(list(globals().keys()) + list(__all__)))
