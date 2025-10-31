"""
Top-level schubmult package initializer.

Dynamically import and re-export all classes defined in all submodules
under this package. Only public classes (names not starting with '_')
that are defined in their module are exported.

This keeps the namespace up-to-date without manually listing every class.
Import errors in individual submodules are ignored (so import-time failures
in a single module won't break importing the whole package).
"""
from __future__ import annotations

import importlib
import inspect
import pkgutil
from typing import List

__all__: List[str] = []

# Optional: keep a static version here if you want
# __version__ = "3.0.2dev1"

# Walk all submodules of this package and import classes
for finder, modname, ispkg in pkgutil.walk_packages(__path__, prefix=__name__ + "."):
    if modname != "scripts":
        mod = importlib.import_module(modname)

    for attr_name, attr_val in vars(mod).items():
        # export only public classes defined in that module
        if attr_name.startswith("_"):
            continue
        if inspect.isclass(attr_val) and getattr(attr_val, "__module__", None) == mod.__name__:
            globals()[attr_name] = attr_val
            __all__.append(attr_name)

# finalize __all__
__all__ = tuple(sorted(__all__))
