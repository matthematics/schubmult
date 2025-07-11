from schubmult.symbolic import DefaultPrinting
from schubmult.utils.perm_utils import add_perm_dict

from .free_algebra import FreeAlgebraElement


def _tup_concat(tup0, tup1, spot):
    if len(tup0) < spot:
        return (*tup0, *tup1)
    return (*tup0[:spot], *tup1, *tup0[spot:])


class SplitFreeAlgebraModule(DefaultPrinting):
    def __init__(self, ring, splitting_spot):
        self._spot = splitting_spot
        self._ring = ring
        self.dtype = type("SplitFreeAlgebraModuleElement", (SplitFreeAlgebraModuleElement,), {"module": self})

    @property
    def zero(self):
        return self.dtype()

    def from_dict(self, dct):
        return self.dtype(dct)

    def new(self, x):
        if isinstance(x, tuple) or isinstance(x, list):
            return self.dtype({tuple(x): 1})
        return self.from_dict(self._ring.new(x))

    def add(self, elem, other):
        return self.from_dict(add_perm_dict(elem, other))

    def lmul(self, other, elem):
        new_elem = self.dtype()
        if isinstance(other, FreeAlgebraElement):
            for tup0, v0 in other.items():
                for tup1, v1 in elem.items():
                    new_elem += self.from_dict({_tup_concat(tup0, tup1, self._spot): v0 * v1})
        return new_elem


class SplitFreeAlgebraModuleElement(DefaultPrinting, dict):
    def __mul__(self, other):
        return other.__rmul__(self)

    def __rmul__(self, other):
        return self.module.lmul(other, self)

    def _sympystr(self, printer):
        return f"M{printer._print(dict(self))}"

    def __add__(self, other):
        if self.module.dtype is other.module.dtype:
            return self.module.add(self, other)
        return NotImplemented
