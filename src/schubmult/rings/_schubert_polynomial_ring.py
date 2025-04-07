from functools import cache

import sympy
from symengine import expand, sympify
from sympy import Add, Basic, Mul
from sympy.core.expr import Expr
from sympy.core.kind import NumberKind
from sympy.printing.str import StrPrinter

import schubmult.rings._utils as utils
import schubmult.schub_lib.double as yz
import schubmult.schub_lib.single as py
from schubmult.perm_lib import (
    Permutation,
    add_perm_dict,
    inv,
)
from schubmult.poly_lib.poly_lib import xreplace_genvars
from schubmult.poly_lib.schub_poly import schubpoly
from schubmult.poly_lib.variables import GeneratingSet, MaskedGeneratingSet
from schubmult.utils.logging import get_logger

## EMULATE POLYTOOLS

_def_printer = StrPrinter({"order": "none"})
# _def_printer = StrPrinter()

logger = get_logger(__name__)

# numpy arrays
# sympy parsing
# quantum

# COPRODUCT


def _varstr(v):
    if v == utils.NoneVar:
        return "NoneVar"
    if v == utils.ZeroVar:
        return "0"
    return f"'{v}'"


# def self._from_dict(coeff_dict):
#     return DoubleSchubertAlgebraElement(coeff_dict)


@cache
def cached_product(u, v, va, vb):
    return {(k, va): xreplace_genvars(x, utils.poly_ring(va), utils.poly_ring(vb)) for k, x in yz.schubmult_double_pair_generic(u, v).items()}


@cache
def cached_positive_product(u, v, va, vb):
    return {(k, va): xreplace_genvars(x, utils.poly_ring(va), utils.poly_ring(vb)) for k, x in yz.schubmult_generic_partial_posify(u, v).items()}


def _mul_schub_dicts(dict1, dict2, best_effort_positive=True):
    by_var = {}

    none_dict = {}
    for k, v in dict1.items():
        if k[1] == utils.NoneVar or k[1] == utils.ZeroVar:
            none_dict[k[0]] = v
        else:
            if k[1] not in by_var:
                by_var[k[1]] = {}
            by_var[k[1]][k[0]] = v

    results = {}
    # import sys

    for _vstr, _dict in by_var.items():
        this_dict = {}
        for k, v in dict2.items():
            for kd, vd in _dict.items():
                did_positive = True
                if best_effort_positive:
                    try:
                        this_dict = add_perm_dict(this_dict, {k1: v1 * v * vd for k1, v1 in cached_positive_product(kd, k[0], _vstr, k[1]).items()})
                    except Exception:
                        logger.debug("Failed to compute")
                        did_positive = False
                    if not did_positive:
                        this_dict = add_perm_dict(this_dict, {k1: v1 * v * vd for k1, v1 in cached_product(kd, k[0], _vstr, k[1]).items()})
        results = add_perm_dict(results, this_dict)

    by_var2 = {}
    none_dict2 = {}
    for k, v in dict2.items():
        if k[1] == utils.NoneVar or k[1] == utils.ZeroVar:
            none_dict2[k[0]] = v
        else:
            if k[1] not in by_var2:
                by_var2[k[1]] = {}
            by_var2[k[1]][k[0]] = v

    for _vstr, _dict in by_var2.items():
        this_dict = {}
        for k, v in none_dict.items():
            if not best_effort_positive:
                this_dict = add_perm_dict(this_dict, {(k1, _vstr): v1 * v for k1, v1 in yz.schubmult_double(_dict, k, utils.poly_ring(_vstr), utils.poly_ring(utils.NoneVar)).items()})
            else:
                this_dict = add_perm_dict(this_dict, {(k1, _vstr): expand(v1) * v for k1, v1 in yz.schubmult_double(_dict, k, utils.poly_ring(_vstr), utils.poly_ring(utils.NoneVar)).items()})
        results = add_perm_dict(results, this_dict)

    none_dict, none_dict2 = sorted([none_dict, none_dict2], key=lambda x: -len(x.keys()))
    for k, v in none_dict2.items():
        results = add_perm_dict(results, {(k1, utils.NoneVar): v1 * v for k1, v1 in py.schubmult_py(none_dict, k).items()})

    return results


class DoubleSchubertAlgebraElement(Expr):
    """Algebra with sympy coefficients
    and a dict basis
    """

    _op_priority = 1e200
    # __slots__ = ("_dict", "_parent")
    _kind = NumberKind
    is_commutative = True
    # is_polynomial = True

    # default_coeff_var = "y"

    def __new__(cls, _dict, genset):
        return DoubleSchubertAlgebraElement.__xnew_cached__(cls, sympy.Dict({k: v for k, v in _dict.items() if expand(v) != 0}), genset)

    @staticmethod
    def __xnew__(_class, _dict, genset):
        obj = Expr.__new__(_class, _dict, genset)
        return obj

    @staticmethod
    @cache
    def __xnew_cached__(_class, _dict, *args, **kwargs):
        return DoubleSchubertAlgebraElement.__xnew__(_class, _dict, *args, **kwargs)

    @property
    def coeff_dict(self):
        return self.args[0]

    @property
    def genset(self):
        return self.args[1]

    def __hash__(self):
        return hash(tuple(self.args))

    # @property
    # def _add_handler(self):
    #     return SchubAdd

    # @property
    # def _mul_handler(self):
    #     return SchubMul

    # confers existing generating set
    def _from_dict(self, _dict):
        return DoubleSchubertAlgebraElement(_dict, self.genset)

    # def _eval_Eq(self, other):
    #     # this will prevent sympy from acting like an idiot
    #     return self.__eq__(other)

    # def _eval_subs(self, old, new):
    #     b_old = sympify(old)
    #     b_new = sympify(new)
    #     result = {}
    #     stuff_to_do = False
    #     lots_of_stuff_to_do = False
    #     if b_new in utils.poly_ring(self._base_var):
    #         stuff_to_do = True
    #     if b_old in utils.poly_ring(self._base_var):
    #         lots_of_stuff_to_do = True
    #     for k, v in self.coeff_dict.items():
    #         if lots_of_stuff_to_do:
    #             poley = sympify(self._from_dict({k: 1}).change_vars(0).expand() * v)
    #             if b_old in poley.free_symbols:
    #                 poley = poley.subs(b_old, b_new)
    #                 new_dict = yz.mult_poly_double({(1, 2): 1}, poley, utils.poly_ring(self._base_var), utils.poly_ring(k[1]))
    #                 new_p = {(koifle, k[1]): voifle for koifle, voifle in new_dict.items()}
    #                 result = add_perm_dict(result, new_p)
    #         elif stuff_to_do:
    #             this_p = self._from_dict({k: v}).change_vars(0)
    #             for kkk, vvv in this_p.coeff_dict.items():
    #                 vvvv = sympify(vvv).subs(b_old, b_new)
    #                 if b_new in sympify(vvvv).free_symbols:
    #                     s_dict = {kkk[0]: 1}
    #                     r_dict = py.mult_poly_py(s_dict, vvvv, utils.poly_ring(self._base_var))
    #                 else:
    #                     r_dict = {kkk[0]: vvvv}
    #                 r_dict = {(kk, 0): voif for kk, voif in r_dict.items()}
    #                 new_p = self._from_dict(r_dict).change_vars(k[1])
    #                 result = add_perm_dict(result, new_p.coeff_dict)
    #         else:
    #             result[k] = result.get(k, 0) + sympify(v).subs(b_old, b_new)
    #     return self._from_dict(result)

    @cache
    def _cached_sympystr(self, printer):
        return printer._print_Add(
            sympy.Add(*[sympy.Mul(self.coeff_dict[k], DSchubPoly(k, self.genset)) for k in sorted(self.coeff_dict.keys(), key=lambda bob: (inv(bob[0]), str(bob[1]), *bob[0]))], evaluate=False),
        )

    def _sympystr(self, printer):
        return self._cached_sympystr(printer)

    # def _eval_simplify(self, *args, measure, **kwargs):
    #     return self._from_dict({k: sympify(sympy.simplify(v, *args, measure=measure, **kwargs)) for k, v in self.coeff_dict.items()})

    def __add__(self, other):
        other = DSx(other)
        return self._from_dict(add_perm_dict(self.coeff_dict, other.coeff_dict))

    def __radd__(self, other):
        other = DSx(other)
        return self._from_dict(add_perm_dict(other.coeff_dict, self.coeff_dict))

    def __sub__(self, other):
        other = DSx(other)
        double_dict = add_perm_dict(self.coeff_dict, {k: -v for k, v in other.coeff_dict.items()})
        return self._from_dict(double_dict)

    def __rsub__(self, other):
        other = DSx(other)
        double_dict = add_perm_dict(other.coeff_dict, {k: -v for k, v in self.coeff_dict.items()})
        return self._from_dict(double_dict)

    def __neg__(self):
        elem = self
        if self.is_Add or self.is_Mul:
            elem = self.doit()
        double_dict = {k: -sympify(v) for k, v in elem.coeff_dict.items()}
        return self._from_dict(double_dict)

    def __mul__(self, other):
        other = DSx(other)
        return self._from_dict(_mul_schub_dicts(self.coeff_dict, other.coeff_dict))

    def __rmul__(self, other):
        other = DSx(other)
        return self._from_dict(_mul_schub_dicts(other.coeff_dict, self.coeff_dict))

    # def equals(self, other):
    #     return self.__eq__(other)

    # def test_equality(self, other, disp=False):
    #     elem1 = self
    #     elem2 = other
    #     done = set()
    #     import sys

    #     for k, v in elem1.coeff_dict.items():
    #         done.add(k)
    #         if expand(v - elem2.coeff_dict.get(k, 0)) != 0:
    #             if disp:
    #                 print(f"{k=} {v=} {elem2.coeff_dict.get(k, 0)=} {expand(v - elem2.coeff_dict.get(k, 0))=}", file=sys.stderr)
    #             return False
    #     for k, v in elem2.coeff_dict.items():
    #         if k in done:
    #             continue
    #         if expand(v - elem1.coeff_dict.get(k, 0)) != 0:
    #             if disp:
    #                 print(f"{k=} {v=} {expand(v - elem1.coeff_dict.get(k, 0))=}", file=sys.stderr)
    #             return False
    #     return True

    # def __eq__(self, other):
    #     if self.is_Add or self.is_Mul:
    #         return self.doit().equals(other)
    #     cv = "y"
    #     elem1 = self
    #     elem2 = other

    #     if not elem1.test_equality(elem2):
    #         elem1_o = elem1.change_vars(cv)
    #         elem2_o = elem2.change_vars(cv)
    #         return elem1_o.test_equality(elem2_o)
    #     return True
    # assert all([k[1] == cv for k in elem1.coeff_dict.keys()])
    # assert all([k[1] == cv for k in elem2.coeff_dict.keys()])

    # def __str__(self):
    #     pieces = []
    #     keys = list(self.coeff_dict.keys())
    #     for k in sorted(keys, key=lambda b: (inv(b[0]), b[1], *b[0])):
    #         v = self.coeff_dict[k]
    #         dvar = "D"
    #         if sympy.expand(v) != 0:
    #             pieces += [
    #                 sympy.Mul(
    #                     v,
    #                     sympy.Symbol(
    #                         f"{dvar}S{DSx._base_var}({list(k[0])}, {_varstr(k[1])})",
    #                         commutative=False,
    #                     )
    #                     if k[0] != Permutation([])
    #                     else 1,
    #                 ),
    #             ]
    #     return sympy.sstr(sympy.Add(*pieces, evaluate=False), order="none")

    # def __repr__(self):
    #     return str(self)
    # def _as_ordered_terms(self, *args, **kwargs):

    @cache
    def change_vars(self, cv):
        result = {}
        for k, v in self.coeff_dict.items():
            result = add_perm_dict(result, {k1: v1 * v for k1, v1 in cached_positive_product(Permutation([]), k[0], cv, k[1]).items()})
        return self._from_dict(result)

    def as_coefficients_dict(self):
        return self.coeff_dict

    # TODO: Masked generating set labels
    def coproduct(self, indices, coeff_var="y", gname1=None, gname2=None):
        result_dict = {}
        if gname1 is None:
            gname1 = f"{self.genset.label}_A"
        if gname2 is None:
            gname2 = f"{self.genset.label}_B"
        gens2 = MaskedGeneratingSet(self.genset, indices)
        gens1 = gens2.complement()
        print(f"{gens1.index_mask=}")
        print(f"{gens2.index_mask=}")
        gens1.set_label(gname1)
        gens2.set_label(gname2)
        for k, v in self.coeff_dict.items():
            key = k[0]
            var_str = k[1]
            # print(f"{var_str=}")
            # print(f"{coeff_var=}")
            coprod_dict = yz.schub_coprod_double(key, indices, utils.poly_ring(var_str), utils.poly_ring(coeff_var))
            # print(f"{coprod_dict=}")
            result_dict = add_perm_dict(result_dict, {((k1, var_str), (k2, coeff_var)): v for (k1, k2), v in coprod_dict.items()})
        result = sympy.Integer(0)
        for ktuple, v in result_dict.items():
            result += sympy.Mul(v, DSchubPoly(ktuple[0], gens1), DSchubPoly(ktuple[1], gens2), evaluate=False)
        return result
        # # will not allow zeros

        # return {k: v for k, v in self.coeff_dict.items() if expand(v) != 0}

    # def normalize_coefficients(self, coeff_var):
    #     return DSx([1, 2], coeff_var) * self

    def expand(self, *_a, **_):
        if isinstance(self, SchubAdd):
            return self.doit().expand()
        if isinstance(self, SchubMul):
            return self.doit().expand()
        return expand(Add(*[v * schubpoly(k[0], self.genset, utils.poly_ring(k[1])) for k, v in self.coeff_dict.items()]))

    def as_polynomial(self):
        return sympy.sympify(Add(*[v * schubpoly(k[0], self.genset, utils.poly_ring(k[1])) for k, v in self.coeff_dict.items()]))


# Atomic Schubert polynomial
class DSchubPoly(DoubleSchubertAlgebraElement):
    is_Atom = True

    def __new__(cls, k, genset):
        return DSchubPoly.__xnew_cached__(cls, k, genset)

    @staticmethod
    def __xnew__(_class, k, genset):
        obj = DoubleSchubertAlgebraElement.__new__(_class, sympy.Dict({(Permutation(k[0]), k[1]): 1}), genset)
        obj._key = k
        return obj

    @staticmethod
    @cache
    def __xnew_cached__(_class, k, genset):
        return DSchubPoly.__xnew__(_class, k, genset)

    def _sympystr(self, printer):
        if self._key[1] == 0 or self._key[1] == utils.NoneVar:
            return printer.doprint(f"S{self.genset.label}({list(self._key[0])})")
        return printer.doprint(f"DS{self.genset.label}({list(self._key[0])}, {_varstr(self._key[1])})")


# None is faster to store
class DoubleSchubertAlgebraElement_basis(Basic):
    def __new__(cls, genset):
        return Basic.__new__(cls, genset)

    @property
    def genset(self):
        return self.args[0]

    def __call__(self, x, cv=None):
        if isinstance(x, list) or isinstance(x, tuple):
            if cv is None:
                cv = "y"
            elem = DoubleSchubertAlgebraElement({(Permutation(x), cv): 1}, self.genset)
        elif isinstance(x, Permutation):
            if cv is None:
                cv = "y"
            elem = DoubleSchubertAlgebraElement({(x, cv): 1}, self.genset)

        elif isinstance(x, DoubleSchubertAlgebraElement):
            if x.is_Add or x.is_Mul:
                return x
            if x.genset == self.genset:
                elem = DoubleSchubertAlgebraElement(x.coeff_dict, self.genset)  # , self)
            else:
                return self(x.expand(), cv)
        else:
            x = sympify(x)
            if cv is None or cv == utils.NoneVar:
                cv = utils.NoneVar
                result = py.mult_poly_py({Permutation([]): 1}, x, self.genset)
            else:
                result = yz.mult_poly_double({Permutation([]): 1}, x, self.genset, utils.poly_ring(cv))
            elem = DoubleSchubertAlgebraElement({(k, cv): v for k, v in result.items()}, self.genset)
        return elem


# def _do_schub_mul(a, b):
#     A = DSx(a)
#     B = DSx(b)
#     return self._from_dict(_mul_schub_dicts(A.coeff_dict, B.coeff_dict))


# def _do_schub_add(a, b):
#     A = DSx(a)
#     B = DSx(b)
#     return self._from_dict(add_perm_dict(A.coeff_dict, B.coeff_dict))


def get_postprocessor(cls):
    if cls is Mul:
        return lambda expr: SchubMul(*expr.args, evaluate=False)  # .doit()
    if cls is Add:
        return lambda expr: SchubAdd(*expr.args, evaluate=False)  # .doit()
    return None


# Basic._constructor_postprocessor_mapping[DoubleSchubertAlgebraElement] = {
#     "Mul": [get_postprocessor(Mul)],
#     "Add": [get_postprocessor(Add)],
# }

# add.register_handlerclass((Expr, SchubAdd), SchubAdd)
# mul.register_handlerclass((Expr, SchubMul), SchubMul)


DSx = DoubleSchubertAlgebraElement_basis(GeneratingSet("x"))


def Sx(x):
    return DSx(x, utils.NoneVar)


DoubleSchubertPolynomial = DoubleSchubertAlgebraElement


class SchubAdd(DoubleSchubertAlgebraElement, Add):
    is_Add = True

    def __new__(cls, *args, evaluate=True, _sympify=True):
        obj = Add.__new__(cls, *args, evaluate=evaluate, _sympify=_sympify)
        if evaluate:
            return obj.doit()
        return obj

    # def doit(self):
    #     ret = self.args[0]
    #     for arg in self.args[1:]:
    #         if arg.is_Add or arg.is_Mul:
    #             arg = arg.doit()
    #         ret = _do_schub_add(ret, arg)
    #     return ret

    # def _sympystr(self, printer):
    #     return _def_printer._print(f"SchubAdd({self.args}")


class SchubMul(DoubleSchubertAlgebraElement, Mul):
    is_Mul = True

    def __new__(cls, *args, evaluate=True, _sympify=True):
        if len(args) == 0:
            return 1
        # args, a, b = Mul.flatten(list(args))
        # if len(args) == 0:
        #     return 1
        obj = Mul.__new__(cls, *args, evaluate=evaluate, _sympify=_sympify)

        if evaluate:
            return obj.doit()
        return obj

    # def doit(self):
    #     ret = self.args[0]
    #     for arg in self.args[1:]:
    #         if arg.is_Add or arg.is_Mul:
    #             arg = arg.doit()
    #         ret = _do_schub_mul(ret, arg)
    #     return ret


# Basic._constructor_postprocessor_mapping[DoubleSchubertAlgebraElement] = {
#     "Mul": [get_postprocessor(Mul)],
#     "Add": [get_postprocessor(Add)],
# }
