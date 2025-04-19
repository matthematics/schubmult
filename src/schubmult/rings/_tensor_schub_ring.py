from functools import cache

import symengine
from symengine import expand
from sympy import Add, Basic, Dict, Expr, Mul, sympify

import schubmult.perm_lib as pl


def _tensor_product_of_dicts(d1, d2):
    ret_dict = {}
    for k1, v1 in d1.items():
        this_dict = {}
        for k2, v2 in d2.items():
            this_dict[(k1, k2)] = v1 * v2
        ret_dict = pl.add_perm_dict(ret_dict, this_dict)
    return Dict(ret_dict)


class TensorBasisElement(Expr):
    is_commutative = False

    def __new__(cls, k1, k2, basis, tensor_symbol=" # "):
        obj = Expr.__new__(cls, k1, k2, basis)
        obj._elem1 = basis.basis1._from_dict({k1: 1})
        obj._elem2 = basis.basis2._from_dict({k2: 1})
        obj._tensor_symbol = tensor_symbol
        return obj

    def _sympystr(self, printer):
        return f"{printer.doprint(self._elem1)}{self._tensor_symbol}{printer.doprint(self._elem2)}"


class TensorAlgebraElement(Expr):
    # tensor ring
    def __new__(cls, _dict, basis):
        return TensorAlgebraElement.__xnew_cached__(cls, Dict(_dict), basis)

    @staticmethod
    def __xnew__(_class, _dict, basis):
        _dict = Dict({k: v for k, v in _dict.items() if expand(v) != 0})
        return Expr.__new__(_class, _dict, basis)

    @property
    def coeff_dict(self):
        return self.args[0]

    @staticmethod
    @cache
    def __xnew_cached__(_class, _dict, basis):
        return TensorAlgebraElement.__xnew__(_class, _dict, basis)

    @property
    def basis(self):
        return self.args[1]

    def __mul__(self, other):
        ret_dict = {}
        for k1, v1 in self.coeff_dict.items():
            for k2, v2 in other.coeff_dict.items():
                dict1 = self.basis.basis1._from_dict({k1[0]: v1 * v2}) * self.basis.basis1._from_dict({k2[0]: 1})
                dict2 = self.basis.basis2._from_dict({k1[1]: 1}) * self.basis.basis2._from_dict({k2[1]: 1})
                ret_dict = pl.add_perm_dict(ret_dict, _tensor_product_of_dicts(dict1.coeff_dict, dict2.coeff_dict))
        return self.basis._from_dict(ret_dict)

    def __add__(self, other):
        return self.basis._from_dict(pl.add_perm_dict(self.coeff_dict, other.coeff_dict))

    def __sub__(self, other):
        return self.basis._from_dict(pl.add_perm_dict(self.coeff_dict, {k: -v for k, v in other.coeff_dict.items()}))

    @cache
    def cached_sympystr(self, printer):
        ret_list = [Mul(v, TensorBasisElement(k[0], k[1], self.basis)) for k, v in self.coeff_dict.items()]
        return printer.doprint(Add(*ret_list))

    def _sympystr(self, printer):
        return self.cached_sympystr(printer)

    def expand(self, **_):
        return sympify(
            symengine.Add(
                *[v * symengine.sympify(self.basis.basis1._from_dict({k[0]: 1}).expand()) * symengine.sympify(self.basis.basis2._from_dict({k[1]: 1}).expand()) for k, v in self.coeff_dict.items()],
            ),
        )


class TensorAlgebraBasis(Basic):
    # tensor ring
    def __new__(cls, basis1, basis2):
        return TensorAlgebraBasis.__xnew_cached__(cls, basis1, basis2)

    @staticmethod
    def __xnew__(_class, basis1, basis2):
        return Basic.__new__(_class, basis1, basis2)

    @property
    def basis1(self):
        return self.args[0]

    @property
    def basis2(self):
        return self.args[1]

    @property
    def coeff_dict(self):
        return self.args[0]

    def _from_dict(self, _dict):
        return TensorAlgebraElement(_dict, self)

    @staticmethod
    @cache
    def __xnew_cached__(_class, basis1, basis2):
        return TensorAlgebraBasis.__xnew__(_class, basis1, basis2)

    def call2(self, *args):
        def calla(*a):
            return TensorAlgebraElement._from_dict(_tensor_product_of_dicts(self.basis1(*args).coeff_dict, self.basis2(*a).coeff_dict))

        return calla

    def __call__(self, *args):
        return self.call2(args)


# def TensorAlgebra_basis(Basic):
#     pass
