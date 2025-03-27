import sympy
from symengine import expand, sympify
from symengine.lib.symengine_wrapper import Expr

import schubmult.schubmult_double as yz
import schubmult.schubmult_py as py
from schubmult.perm_lib import add_perm_dict, permtrim

from ._utils import poly_ring

# numpy arrays
# sympy parse


def _mul_schub_dicts(dict1, dict2):
    results = {}

    for k, v in dict2.items():
        results = add_perm_dict(results, {k: v0 * v for k, v0 in py.schubmult(dict1, k).items()})
    return results


class DictAlgebraElement:
    """Algebra with sympy coefficients
    and a dict basis
    """

    def __init__(self, _dict, parent):
        self._dict = _dict
        self._parent = parent

    def __add__(self, other):
        try:
            return DictAlgebraElement(add_perm_dict(self._dict, self._parent(other)._dict), self._parent)
        except Exception:
            return NotImplemented

    def __radd__(self, other):
        try:
            return DictAlgebraElement(add_perm_dict(self._parent(other)._dict, self._dict), self._parent)
        except Exception:
            return NotImplemented

    def __sub__(self, other):
        try:
            return DictAlgebraElement(add_perm_dict(self._dict, {k: -v for k, v in self._parent(other)._dict.items()}), self._parent)
        except Exception:
            return NotImplemented

    def __rsub__(self, other):
        try:
            return DictAlgebraElement(add_perm_dict(self._parent(other)._dict, {k: -v for k, v in self._dict.items()}), self._parent)
        except Exception:
            return NotImplemented

    def __neg__(self):
        try:
            return DictAlgebraElement({k: -v for k, v in self._dict.items()}, self._parent)
        except Exception:
            return NotImplemented

    def __mul__(self, other):
        try:
            return DictAlgebraElement(_mul_schub_dicts(self._dict, self._parent(other)._dict), self._parent)
        except Exception:
            return NotImplemented

    def __rmul__(self, other):
        try:
            return DictAlgebraElement(_mul_schub_dicts(self._parent(other)._dict, self._dict), self._parent)
        except Exception:
            return NotImplemented

    def __eq__(self, other):
        elem1 = self
        elem2 = self._parent(other)
        done = set()
        for k, v in elem1._dict.items():
            done.add(k)
            if expand(v - elem2._dict.get(k, 0)) != 0:
                return False
        for k, v in elem2._dict.items():
            if k in done:
                continue
            if expand(v - elem1._dict.get(k, 0)) != 0:
                return False
        return True

    def __str__(self):
        pieces = []
        keys = list(self._dict.keys())
        for k in sorted(keys):
            v = self._dict[k]
            if expand(v) != 0:
                pieces += [sympy.Mul(v, sympy.Symbol(f"S{self._parent._base_var}({list(k[0])})", commutative=False))]
        return sympy.sstr(sympy.Add(*pieces), order=lambda order: pieces)  # use sstr

    def __repr__(self):
        return self.__str__()

    def __expand__(self):
        return self.expand()

    # def _sympify(self):
    #     return self

    def expand(self, deep=True):  # noqa: ARG002
        ret = 0
        for k, v in self._dict.items():
            ret += yz.schubmult({(1, 2): v}, k, poly_ring(self._parent._base_var), [0 for i in range(100)]).get((1, 2), 0)
        return ret

    def as_coefficients_dict(self):
        # will not allow zeros
        return {k: v for k, v in self._dict.items() if expand(v) != 0}


class DictAlgebraElement_basis:
    def __init__(self, base_var="x"):
        self._base_var = base_var

    def __call__(self, x):
        if isinstance(x, list) or isinstance(x, tuple):
            elem = DictAlgebraElement({tuple(permtrim(list(x))): 1}, self)
        elif isinstance(x, DictAlgebraElement):
            if x._parent._base_var == self._base_var:
                elem = DictAlgebraElement(x._dict, self)
            else:
                return self(x.expand())
        else:
            try:
                result = py.mult_poly({(1, 2): 1}, sympify(x), poly_ring(self._base_var))  # this will raise an error if no good
                elem = DictAlgebraElement(result, self)
            except Exception:
                # return NotImplemented
                return NotImplemented
        return elem


Schub = DictAlgebraElement_basis()
SchubertPolynomial = DictAlgebraElement
