import sympy
from symengine import Expr, expand, sympify

import schubmult.rings._schubert_polynomial_ring as spr
import schubmult.schubmult_double as yz
import schubmult.schubmult_py as py
import symengine.lib.symengine_wrapper as syme_wrap
from schubmult.perm_lib import add_perm_dict, permtrim

from ._utils import poly_ring

# numpy arrays
# sympy parsing
# quantum


def _mul_schub_dicts(dict1, dict2):
    by_var = {}

    for k, v in dict1.items():
        if k[1] not in by_var:
            by_var[k[1]] = {}
        by_var[k[1]][k[0]] = v

    results = {}

    for _vstr, _dict in by_var.items():
        this_dict = {}
        for k, v in dict2.items():
            if k[1] is None and _vstr is None:
                this_dict = add_perm_dict(this_dict, {(k1, None): v1 * v for k1, v1 in py.schubmult(_dict, k[0])})
            else:
                this_dict = add_perm_dict(this_dict, {(k1, _vstr): v1 * v for k1, v1 in yz.schubmult(_dict, k[0], poly_ring(_vstr), poly_ring(k[1])).items()})
        results.update(this_dict)
    return results

from functools import cache
class DoubleDictAlgebraElement(Expr):
    """Algebra with sympy coefficients
    and a dict basis
    """

    _op_priority = 1e200
    __slots__ = ("_dict", "_parent")

    # def __hash__(self):
    #     return hash(sympy.Dict(self._dict))

    def __new__(cls, _dict, parent):
        print("Prifbonk")
        return DoubleDictAlgebraElement.__xnew_cached__(sympy.Dict(_dict), parent)

    @classmethod
    def __xnew__(cls, _dict, parent):
        obj = Expr.__new__(cls)
        obj._dict = _dict
        obj._parent = parent
        print("ASFBIASF")
        return obj

    @classmethod
    @cache
    def __xnew_cached__(cls, _dict, parent):
        print(f"{type(_dict)=} {hash(parent)}")
        return DoubleDictAlgebraElement.__xnew__(_dict, parent)

    @classmethod
    @cache
    def __xnew_sympy_cached__(cls, _dict, parent):
        print("Plefflenoff gaboopa")
        return DoubleDictAlgebraElement.__xnew_sympy__(_dict, parent)

    def _eval_simplify_(self):
        print(f"{self=} preffdoffle")
        return self

    @classmethod
    def __xnew_sympy__(cls, _dict, parent):
        obj = syme_wrap.Expr.__new__(cls)
        obj._dict = _dict
        obj._parent = parent
        print("ASFBIAS brofool F")
        return obj

    def _symengine_(self):
        print("piffwaffle")
        return DoubleDictAlgebraElement.__xnew_cached__(self._dict, self._parent)
    
    # def _sympy_(self):
    #     print("Profdonket")
    #     # print(f"{self=}")
    #     return DoubleDictAlgebraElement.__xnew_sympy_cached__(self._dict, self._parent)

    # def _symplify_(self):
    #     print("Frababa")
    # def __init__(self, _dict, parent):
    #     self._dict = _dict
    #     self._parent = parent

    is_Atom = True

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
        # print("FERPE")
        one = self._parent([1, 2])
        elem1 = one * self
        elem2 = one * self._parent(other)
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
        # print("ASFJAdsajdSJ")

        pieces = []
        keys = list(self._dict.keys())
        for k in sorted(keys):
            v = self._dict[k]
            if expand(v) != 0:
                pieces += [
                    sympy.Mul(
                        v,
                        sympy.Symbol(
                            ((f"S{self._parent._base_var}" if k[1] is None else f"DS{self._parent._base_var}") + (f"({list(k[0])})" if k[1] is None else f"({list(k[0])}, '{k[1]}')")),
                            commutative=False,
                        ),
                    ),
                ]
        return sympy.sstr(sympy.Add(*pieces), order=lambda order: pieces)  # use sstr

    def __repr__(self):
        # print("ASFJAdsajdSJ")

        return self.__str__()

    def as_coefficients_dict(self):
        # will not allow zeros
        # print("ASFJAdsajdSJ")

        return {k: v for k, v in self._dict.items() if expand(v) != 0}

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
            ret = DoubleDictAlgebraElement({k: expand(v) for k, v in self._dict.items()}, self._parent)
        return ret


class DoubleDictAlgebraElement_basis:
    def __init__(self, base_var="x", coeff_var=None):
        # self._dict = _dict
        self._base_var = base_var
        self._coeff_var = coeff_var if coeff_var else "y"

    def __call__(self, *x):
        # print(f"Profilingui {self=} {x=}")
        cv = ""
        if len(x) == 1:
            x = x[0]
        elif len(x) == 2 and (isinstance(x[0], list) or isinstance(x[0], tuple)):
            x, cv = x[0], x[1]
        if cv == "":
            cv = self._coeff_var
        if isinstance(x, list) or isinstance(x, tuple):
            elem = DoubleDictAlgebraElement({(tuple(permtrim(list(x))), cv): 1}, self)
        elif isinstance(x, DoubleDictAlgebraElement):
            if x._parent._base_var == self._base_var:
                elem = DoubleDictAlgebraElement(x._dict, self)
            else:
                return self(x.expand())
        elif isinstance(x, spr.SchubertPolynomial):
            if x._parent._base_var == self._base_var:
                elem_dict = {(k, None): v for k, v in x._dict.items()}
                elem = DoubleDictAlgebraElement(elem_dict, self)
            else:
                return self(x.expand())
        else:
            result = yz.mult_poly({(1, 2): 1}, sympify(x), poly_ring(self._base_var), poly_ring(self._coeff_var))
            elem = DoubleDictAlgebraElement({(k, self._coeff_var): v for k, v in result.items()}, self)
        return elem

    # def _coerce_map_from(self, S):
    #     if isinstance(S, type(DSchub)):
    #         return True
    #     if isinstance(S, type(Schub)):
    #         return True
    #     if isinstance(S, Expr):
    #         return True
    #     return False


DSchub = DoubleDictAlgebraElement_basis()
DoubleSchubertPolynomial = DoubleDictAlgebraElement
