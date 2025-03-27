from functools import cache

import sympy

import schubmult.rings._schubert_polynomial_ring as spr
import schubmult.schubmult_double as yz
import schubmult.schubmult_py as py
from schubmult.perm_lib import add_perm_dict, inv, permtrim

from ._utils import NoneVar, ZeroVar, poly_ring

# numpy arrays
# sympy parsing
# quantum
# none is faster
# coeff_var = None
# singleton
# class NoneVarType(sympy.Basic):
#     def __new__(cls, *args):
#         return NoneVarType.__xnew_cached__(cls, args)

#     def __xnew__(cls, *Args):
#         obj = sympy.Basic.__new__(*args)
#         return obj


# def is_nonevar(v):
#     return isinstance(v, NoneVarType)

# _NoneVar = NoneVarType()

# sympy.init_printing(order='none')


def _varstr(v):
    if v == NoneVar:
        return "NoneVar"
    if v == ZeroVar:
        return "0"
    return f"'{v}'"


def _mul_schub_dicts(dict1, dict2):
    by_var = {}

    none_dict = {}
    for k, v in dict1.items():
        if k[1] == NoneVar:
            none_dict[k[0]] = v
        else:
            if k[1] not in by_var:
                by_var[k[1]] = {}
            by_var[k[1]][k[0]] = v

    results = {}

    for _vstr, _dict in by_var.items():
        this_dict = {}
        for k, v in dict2.items():
            print(f"{k=}")
            this_dict = add_perm_dict(this_dict, {(k1, k[1]): v1 * v for k1, v1 in yz.schubmult(_dict, k[0], poly_ring(_vstr), poly_ring(k[1])).items()})
        results.update(this_dict)

    by_var2 = {}
    none_dict2 = {}
    for k, v in dict2.items():
        if k[1] == NoneVar:
            none_dict2[k[0]] = v
        else:
            if k[1] not in by_var2:
                by_var2[k[1]] = {}
            by_var2[k[1]][k[0]] = v

    for _vstr, _dict in by_var2.items():
        this_dict = {}
        for k, v in none_dict.items():
            this_dict = add_perm_dict(this_dict, {(k1, _vstr): v1 * v for k1, v1 in yz.schubmult(_dict, k, poly_ring(_vstr), poly_ring(NoneVar)).items()})
        results = add_perm_dict(results, this_dict)

    this_dict = {}
    for k, v in none_dict2.items():
        this_dict = add_perm_dict(this_dict, {(k1, NoneVar): v1 * v for k1, v1 in py.schubmult(none_dict, k).items()})
    results.update(this_dict)

    return results


class DoubleDictAlgebraElement(sympy.Expr):
    """Algebra with sympy coefficients
    and a dict basis
    """

    _op_priority = 1e200
    __slots__ = ("_dict", "_parent")

    # default_coeff_var = "y"

    def __new__(cls, _dict, parent):
        return DoubleDictAlgebraElement.__xnew_cached__(sympy.Dict(_dict), parent)

    @classmethod
    def __xnew__(cls, _dict, parent):
        obj = sympy.Expr.__new__(cls)
        obj._dict = _dict
        obj._parent = parent
        return obj

    @classmethod
    @cache
    def __xnew_cached__(cls, _dict, parent):
        return DoubleDictAlgebraElement.__xnew__(_dict, parent)

    def _symengine_(self):
        return NotImplemented

    def _eval_simplify_(self):
        return self

    def __add__(self, other):
        # print("ASFJASJ")
        return DoubleDictAlgebraElement(add_perm_dict(self._dict, self._parent(other)._dict), self._parent)

    def __radd__(self, other):
        # print("ASFJAdsajdSJ")

        return DoubleDictAlgebraElement(add_perm_dict(self._parent(other)._dict, self._dict), self._parent)

    def __sub__(self, other):
        # print("ASFJAdsajdSJ")

        return DoubleDictAlgebraElement(add_perm_dict(self._dict, {k: -v for k, v in self._parent(other)._dict.items()}), self._parent)

    def __rsub__(self, other):
        # print("ASFJAdsajdSJ")

        return DoubleDictAlgebraElement(add_perm_dict(self._parent(other)._dict, {k: -v for k, v in self._dict.items()}), self._parent)

    def __neg__(self):
        # print("ASFJAdsajdSJ")

        return DoubleDictAlgebraElement({k: -v for k, v in self._dict.items()}, self._parent)

    def __mul__(self, other):
        # print("ASFJAdsajdSJ")

        return DoubleDictAlgebraElement(_mul_schub_dicts(self._dict, self._parent(other)._dict), self._parent)

    def __rmul__(self, other):
        # print("ASFJAdsajdSJ")

        return DoubleDictAlgebraElement(_mul_schub_dicts(self._parent(other)._dict, self._dict), self._parent)

    def __eq__(self, other):
        elem1 = self._parent(self, "y")  # count vars?
        elem2 = self._parent(other, "y")
        done = set()
        for k, v in elem1._dict.items():
            done.add(k)
            if sympy.expand(v - elem2._dict.get(k, 0)) != 0:
                return False
        for k, v in elem2._dict.items():
            if k in done:
                continue
            if sympy.expand(v - elem1._dict.get(k, 0)) != 0:
                return False
        return True

    def __str__(self):
        pieces = []
        keys = list(self._dict.keys())
        for k in sorted(keys, key=lambda b: (inv(b[0]), b[1], *b[0])):
            v = self._dict[k]
            dvar = "D"
            if sympy.expand(v) != 0:
                pieces += [
                    sympy.Mul(
                        v,
                        sympy.Symbol(
                            f"{dvar}S{self._parent._base_var}({list(k[0])}, {_varstr(k[1])})",
                            commutative=False,
                        )
                        if k[0] != (1, 2)
                        else 1,
                    ),
                ]
        return sympy.sstr(sympy.Add(*pieces, evaluate=False), order="none")

    def __repr__(self):
        return str(self)

    def as_coefficients_dict(self):
        # will not allow zeros
        # print("ASFJAdsajdSJ")

        return {k: v for k, v in self._dict.items() if sympy.expand(v) != 0}

    def normalize_coefficients(self, coeff_var):
        return self._parent([1, 2], coeff_var) * self

    # def schub_coeff(self, perm):
    #     return self._dict.get(tuple(permtrim(perm)), 0)

    def expand(self, deep=True):
        if deep:
            ret = 0
            keys = list(self._dict.keys())
            for k in sorted(keys):
                v = self._dict[k]
                ret += yz.schubmult({(1, 2): v}, k[0], poly_ring(self._parent._base_var), poly_ring(k[1])).get((1, 2), 0)
        else:
            ret = DoubleDictAlgebraElement({k: sympy.expand(v) for k, v in self._dict.items()}, self._parent)
        return sympy.sympify(ret)


# None is faster to store
class DoubleDictAlgebraElement_basis:
    coeff_varname = "y"

    def __init__(self, base_var="x"):
        # self._dict = _dict
        self._base_var = base_var
        # self._coeff_var = coeff_var if coeff_var else "y"

    def __call__(self, x, cv=None):
        if isinstance(x, list) or isinstance(x, tuple):
            if cv is None:
                cv = "y"
            elem = DoubleDictAlgebraElement({(tuple(permtrim(list(x))), cv): 1}, self)
        elif isinstance(x, DoubleDictAlgebraElement):
            if x._parent._base_var == self._base_var:
                elem = DoubleDictAlgebraElement(x._dict, self)
            else:
                return self(x.expand(), cv)
        elif isinstance(x, spr.SchubertPolynomial):
            if x._parent._base_var == self._base_var:
                elem_dict = {(x, NoneVar): v for k, v in x._dict.items()}
                elem = DoubleDictAlgebraElement(elem_dict, self)
                if cv is not None:
                    elem = self([1, 2], cv) * elem
            else:
                return self(x.expand(), cv)
        else:
            if cv is None:
                cv = NoneVar
                result = py.mult_poly({(1, 2): 1}, x, poly_ring(self._base_var))
            else:
                result = yz.mult_poly({(1, 2): 1}, x, poly_ring(self._base_var), poly_ring(cv))
            elem = DoubleDictAlgebraElement({(k, cv): v for k, v in result.items()}, self)
        return elem

    # def _coerce_map_from(self, S):
    #     if isinstance(S, type(DSchub)):
    #         return True
    #     if isinstance(S, type(Schub)):
    #         return True
    #     if isinstance(S, Expr):
    #         return True
    #     return False


DSx = DoubleDictAlgebraElement_basis()
DoubleSchubertPolynomial = DoubleDictAlgebraElement
