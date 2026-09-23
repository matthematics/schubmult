"""`LazyAttr`: a stand-in for a module attribute that imports the module on first real use."""

import importlib
import sys


class LazyAttr:
    """Proxy for ``getattr(import_module(modname), attr)``, resolved on first call, attribute access,
    or use as a base class. ``isinstance``/``issubclass`` checks resolve it only once ``modname`` has
    been imported, since no instance of a class can exist before its module is loaded.
    """

    __slots__ = ("_attr", "_modname", "_obj")

    def __init__(self, modname, attr):
        self._modname = modname
        self._attr = attr
        self._obj = None

    def _resolve(self):
        if self._obj is None:
            self._obj = getattr(importlib.import_module(self._modname), self._attr)
        return self._obj

    def __call__(self, *args, **kwargs):
        return self._resolve()(*args, **kwargs)

    def __getattr__(self, name):
        return getattr(self._resolve(), name)

    def __instancecheck__(self, instance):
        return self._modname in sys.modules and isinstance(instance, self._resolve())

    def __subclasscheck__(self, subclass):
        return self._modname in sys.modules and issubclass(subclass, self._resolve())

    def __mro_entries__(self, bases):
        return (self._resolve(),)

    def __repr__(self):
        return f"<lazy {self._modname}.{self._attr}>"


def lazy_from(modname, *names):
    """``LazyAttr`` proxies for ``from modname import *names``."""
    return tuple(LazyAttr(modname, name) for name in names)
