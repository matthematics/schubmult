# from abc import ABC
from symengine import Add, Mul, Pow, symarray, Basic, expand, Symbol, sympify
from schubmult.schubmult_py import schubmult, mult_poly
from schubmult.perm_lib import add_perm_dict, permtrim
from ._utils import poly_ring

# _def_gens = tuple(symarray("x", 100).tolist())
# _def_cof

class DictAlgebraElement:
    """Algebra with sympy coefficients
    and a dict basis
    """

    def __init__(self, _dict = None,  base_var = "x"):
        self._dict = _dict
        self._base_var = base_var
        # self._base_var = poly_ring(self._base_var)


    def __add__(self, other):
        if isinstance(other, DictAlgebraElement):
            return DictAlgebraElement(add_perm_dict(self._dict, other._dict), self._base_var)
        else:
            return DictAlgebraElement(add_perm_dict(self._dict, self(other)._dict), self._base_var)

    def __sub__(self, other):
        if isinstance(other, DictAlgebraElement):
            return self.__add__(DictAlgebraElement({k: -v for k, v in other._dict.items()}, self._base_var))
        else:
            return self.__add__(DictAlgebraElement({k: -v for k, v in self(other)._dict.items()}, self._base_var))

    def __neg__(self):
        return DictAlgebraElement({k: -v for k, v in self._dict.items()}, self._base_var)

    def _sympify_(self):
        raise TypeError("You cannot sympify this")

    def __mul__(self, other):
        elem = other
        if not isinstance(other, DictAlgebraElement):
            elem = self(other)

        ret = {}
        for k, v in elem._dict.items():
            ret = add_perm_dict(ret, schubmult({k0: v0 * v for k0, v0 in self._dict.items()}, k))
        return DictAlgebraElement(ret, self._base_var)

    def __str__(self):
        pieces = []
        for k, v in self._dict.items():
            if expand(v) != 0:
                pieces += [v * Symbol(f"S{self._base_var}({list(k)})")]
        return str(Add(*pieces))

    def __repr__(self):
        return self.__str__()


class DictAlgebraElement_basis:
    def __init__(self, base_var = "x"):
        self._base_var = base_var

    def __call__(self, x):
        if isinstance(x, list) or isinstance(x, tuple):
            # checking the input to avoid symmetrica crashing Sage, see trac 12924
            elem = DictAlgebraElement({tuple(permtrim(list(x))): 1}, tuple(symarray("x", 100).tolist()))
        else:
            result = mult_poly(
                {(1, 2): 1},
                sympify(x),
                poly_ring(self._base_var),
            )
            elem = DictAlgebraElement(result, tuple(symarray("x", 100).tolist()))
        return elem

Schub = DictAlgebraElement_basis()
SchubertPolynomial = DictAlgebraElement
