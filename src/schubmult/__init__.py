import inspect
import pkgutil
from typing import List

__version__ = "4.0.0dev"

"""
Top-level schubmult package initializer.

Dynamically import and re-export all classes defined in all submodules
under this package. Only public classes (names not starting with '_')
that are defined in their module are exported.

This variant avoids importlib and uses the builtin __import__ to load modules.
"""

__all__: List[str] = []

# Walk all submodules of this package and import classes
for finder, modname, ispkg in pkgutil.walk_packages(__path__, prefix=__name__ + "."):
    # skip any module under the schubmult.scripts package
    if modname.startswith(__name__ + ".scripts"):
        continue

    try:
        # use builtin __import__ (no importlib) and request the module object
        mod = __import__(modname, fromlist=["*"])
    except Exception:
        # ignore modules that fail to import so package import remains robust
        continue

    for attr_name, attr_val in vars(mod).items():
        # export only public classes defined in that module
        if attr_name.startswith("_"):
            continue
        if inspect.isclass(attr_val) and getattr(attr_val, "__module__", None) == mod.__name__:
            globals()[attr_name] = attr_val
            __all__.append(attr_name)

# finalize __all__
__all__ = tuple(sorted(set(__all__)))
