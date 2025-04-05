# to encourage development

from functools import cache

import sympy
from symengine import expand, sympify
from sympy import Add, Basic
from sympy.core.expr import Expr
from sympy.core.kind import NumberKind
from sympy.printing.str import StrPrinter

import schubmult.rings._utils as utils
import schubmult.schub_lib.quantum as py
import schubmult.schub_lib.quantum_double as yz
from schubmult.perm_lib import (
    Permutation,
    add_perm_dict,
    inv,
)
from schubmult.poly_lib.poly_lib import xreplace_genvars

#from schubmult.rings._schubert_polynomial_ring import DoubleSchubertAlgebraElement
from schubmult.schub_lib.quantum_double import schubpoly_quantum
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


def _from_double_dict(_doubledict):
    return QuantumDoubleSchubertAlgebraElement(_doubledict)


@cache
def cached_product(u, v, va, vb):
    return {(k, va): xreplace_genvars(x, utils.poly_ring(va), utils.poly_ring(vb)) for k, x in yz.schubmult_q_double_pair_generic(u, v).items()}


@cache
def cached_positive_product(u, v, va, vb):
    return {(k, va): xreplace_genvars(x, utils.poly_ring(va), utils.poly_ring(vb)) for k, x in yz.schubmult_q_generic_partial_posify(u, v).items()}


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
                this_dict = add_perm_dict(this_dict, {(k1, _vstr): v1 * v for k1, v1 in yz.schubmult_q_double_fast(_dict, k, utils.poly_ring(_vstr), utils.poly_ring(utils.NoneVar)).items()})
            else:
                this_dict = add_perm_dict(this_dict, {(k1, _vstr): expand(v1) * v for k1, v1 in yz.schubmult_q_double_fast(_dict, k, utils.poly_ring(_vstr), utils.poly_ring(utils.NoneVar)).items()})
        results = add_perm_dict(results, this_dict)

    none_dict, none_dict2 = sorted([none_dict, none_dict2], key=lambda x: -len(x.keys()))
    for k, v in none_dict2.items():
        results = add_perm_dict(results, {(k1, utils.NoneVar): v1 * v for k1, v1 in py.schubmult_q_fast(none_dict, k).items()})

    return results


class QuantumDoubleSchubertAlgebraElement(Expr):
    """Algebra with sympy coefficients
    and a dict basis
    """
    is_quantum = True
    _base_var = "x"

    _op_priority = 1e200
    # __slots__ = ("_dict", "_parent")
    _kind = NumberKind
    is_commutative = True
    # is_polynomial = True

    # default_coeff_var = "y"

    def __new__(cls, _dict, *args, **kwargs):
        # for k in _dict.keys():
        #     if not isinstance(k[0], Permutation):
        #         raise TypeError(f"Key {k[0]} is not a permutation")
        return QuantumDoubleSchubertAlgebraElement.__xnew_cached__(cls, sympy.Dict({k: v for k, v in _dict.items() if expand(v) != 0}), *args, **kwargs)

    def __hash__(self):
        return hash(tuple(self.args))

    # @property
    # def _add_handler(self):
    #     return SchubAdd

    # @property
    # def _mul_handler(self):
    #     return SchubMul

    def _eval_Eq(self, other):
        # this will prevent sympy from acting like an idiot
        return self.__eq__(other)

    def _eval_subs(self, old, new):
        b_old = sympify(old)
        b_new = sympify(new)
        result = {}
        stuff_to_do = False
        lots_of_stuff_to_do = False
        if b_new in utils.poly_ring(self._base_var):
            stuff_to_do = True
        if b_old in utils.poly_ring(self._base_var):
            lots_of_stuff_to_do = True
        for k, v in self._doubledict.items():
            if lots_of_stuff_to_do:
                poley = sympify(_from_double_dict({k: 1}).change_vars(0).expand() * v)
                if b_old in poley.free_symbols:
                    poley = poley.subs(b_old, b_new)
                    new_dict = yz.mult_poly_double({(1, 2): 1}, poley, utils.poly_ring(self._base_var), utils.poly_ring(k[1]))
                    new_p = {(koifle, k[1]): voifle for koifle, voifle in new_dict.items()}
                    result = add_perm_dict(result, new_p)
            elif stuff_to_do:
                this_p = _from_double_dict({k: v}).change_vars(0)
                for kkk, vvv in this_p._doubledict.items():
                    vvvv = sympify(vvv).subs(b_old, b_new)
                    if b_new in sympify(vvvv).free_symbols:
                        s_dict = {kkk[0]: 1}
                        r_dict = py.mult_poly_py(s_dict, vvvv, utils.poly_ring(self._base_var))
                    else:
                        r_dict = {kkk[0]: vvvv}
                    r_dict = {(kk, 0): voif for kk, voif in r_dict.items()}
                    new_p = _from_double_dict(r_dict).change_vars(k[1])
                    result = add_perm_dict(result, new_p._doubledict)
            else:
                result[k] = result.get(k, 0) + sympify(v).subs(b_old, b_new)
        return _from_double_dict(result)

    @staticmethod
    def __xnew__(_class, _dict, *args, **kwargs):
        obj = Expr.__new__(_class, _dict)
        obj._doubledict = _dict
        # obj._print_sum = Add(*[sympy.Mul(_dict[k], DSchubSymbol(QDSx._base_var, k)) for k in sorted(_dict.keys(), key=lambda bob: (inv(bob[0]), str(bob[1]), *bob[0]))], evaluate=False)
        return obj

    @cache
    def _cached_sympystr(self, printer):
        # return _def_printer.doprint(self._print_sum)
        return printer._print_Add(
            sympy.Add(*[sympy.Mul(self._doubledict[k], QDSchubPoly(k)) for k in sorted(self._doubledict.keys(), key=lambda bob: (inv(bob[0]), str(bob[1]), *bob[0]))], evaluate=False)
        )

    def _sympystr(self, printer):
        return self._cached_sympystr(printer)

    @staticmethod
    @cache
    def __xnew_cached__(_class, _dict, *args, **kwargs):
        return QuantumDoubleSchubertAlgebraElement.__xnew__(_class, _dict, *args, **kwargs)

    def _symengine_(self):
        return NotImplemented

    def _eval_simplify(self, *args, measure, **kwargs):
        boible = _from_double_dict({k: sympify(sympy.simplify(v, *args, measure=measure, **kwargs)) for k, v in self._doubledict.items()})
        return boible

    def __add__(self, other):
        return _from_double_dict(add_perm_dict(self._doubledict, QDSx(other)._doubledict))

    def __radd__(self, other):

        return _from_double_dict(add_perm_dict(QDSx(other)._doubledict, self._doubledict))

    def __sub__(self, other):
        double_dict = add_perm_dict(self._doubledict, {k: -v for k, v in QDSx(other)._doubledict.items()})
        return _from_double_dict(double_dict)

    def __rsub__(self, other):
        double_dict = add_perm_dict(QDSx(other)._doubledict, {k: -v for k, v in self._doubledict.items()})
        return _from_double_dict(double_dict)        

    def __neg__(self):
        elem = self
        if self.is_Add or self.is_Mul:
            elem = self.doit()
        double_dict = {k: -sympify(v) for k, v in elem._doubledict.items()}
        return _from_double_dict(double_dict)

    def __mul__(self, other):
        return _from_double_dict(_mul_schub_dicts(self._doubledict, QDSx(other)._doubledict))

    def __rmul__(self, other):
        return _from_double_dict(_mul_schub_dicts(QDSx(other)._doubledict, self._doubledict))

    def equals(self, other):
        return self.__eq__(other)

    def test_equality(self, other, disp=False):
        elem1 = self
        elem2 = other
        done = set()
        import sys

        for k, v in elem1._doubledict.items():
            done.add(k)
            if expand(v - elem2._doubledict.get(k, 0)) != 0:
                if disp:
                    print(f"{k=} {v=} {elem2._doubledict.get(k, 0)=} {expand(v - elem2._doubledict.get(k, 0))=}", file=sys.stderr)
                return False
        for k, v in elem2._doubledict.items():
            if k in done:
                continue
            if expand(v - elem1._doubledict.get(k, 0)) != 0:
                if disp:
                    print(f"{k=} {v=} {expand(v - elem1._doubledict.get(k, 0))=}", file=sys.stderr)
                return False
        return True

    def __eq__(self, other):
        if self.is_Add or self.is_Mul:
            return self.doit().equals(other)
        cv = "y"
        elem1 = self
        elem2 = QDSx(other)

        if not elem1.test_equality(elem2):
            elem1_o = elem1.change_vars(cv)
            elem2_o = elem2.change_vars(cv)
            return elem1_o.test_equality(elem2_o)
        return True
        # assert all([k[1] == cv for k in elem1._doubledict.keys()])
        # assert all([k[1] == cv for k in elem2._doubledict.keys()])

    # def __str__(self):
    #     pieces = []
    #     keys = list(self._doubledict.keys())
    #     for k in sorted(keys, key=lambda b: (inv(b[0]), b[1], *b[0])):
    #         v = self._doubledict[k]
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
        for k, v in self._doubledict.items():
            result = add_perm_dict(result, {k1: v1 * v for k1, v1 in cached_positive_product(Permutation([]), k[0], cv, k[1]).items()})
        return _from_double_dict(result)

    def as_coefficients_dict(self):
        # will not allow zeros

        return {k: v for k, v in self._doubledict.items() if expand(v) != 0}

    def normalize_coefficients(self, coeff_var):
        return QDSx([1, 2], coeff_var) * self

    def expand(self, *_a, **_):
        # if isinstance(self, SchubAdd):
        #     return self.doit().expand()
        # if isinstance(self, SchubMul):
        #     return self.doit().expand()
        return expand(Add(*[v * schubpoly_quantum(k[0], utils.poly_ring(QDSx._base_var), utils.poly_ring(k[1])) for k, v in self._doubledict.items()]))

    def as_polynomial(self):
        return sympy.sympify(Add(*[v * schubpoly_quantum(k[0], utils.poly_ring(QDSx._base_var), utils.poly_ring(k[1])) for k, v in self._doubledict.items()]))


# TODO: not a noncommutative symbol, something else
# Atomic Schubert polynomial
class QDSchubPoly(QuantumDoubleSchubertAlgebraElement):
    is_Atom = True

    def __new__(cls, k, *args, **kwargs):
        return QDSchubPoly.__xnew_cached__(k, *args, **kwargs)

    @classmethod
    def __xnew__(cls, k):
        obj = QuantumDoubleSchubertAlgebraElement.__new__(cls, sympy.Dict({(Permutation(k[0]), k[1]): 1}))
        obj._perm = k[0]
        obj._coeff_var = k[1]
        # obj._base_var = base_var
        return obj

    @classmethod
    @cache
    def __xnew_cached__(cls, k):
        return QDSchubPoly.__xnew__(k)

    def _sympystr(self, printer):
        if self._coeff_var == 0 or self._coeff_var == utils.NoneVar:
            return printer.doprint(f"QSx({list(self._perm)})")
        return printer.doprint(f"QDSx({list(self._perm)}, {_varstr(self._coeff_var)})")


# None is faster to store
class QuantumDoubleSchubertAlgebraElement_basis(Basic):
    coeff_varname = "y"

    def __init__(self, base_var="x"):
        # self._doubledict = _dict
        self._base_var = base_var
        # self._coeff_var = coeff_var if coeff_var else "y"

    def __call__(self, x, cv=None):
        if isinstance(x, list) or isinstance(x, tuple):
            if cv is None:
                cv = "y"
            elem = QuantumDoubleSchubertAlgebraElement({(Permutation(x), cv): 1})
        elif isinstance(x, Permutation):
            if cv is None:
                cv = "y"
            elem = QuantumDoubleSchubertAlgebraElement({(x, cv): 1})
        # elif isinstance(x, spr.SchubertPolynomial):
        #     if x._parent._base_var == self._base_var:
        #         elem_dict = {(x, utils.NoneVar): v for k, v in x._doubledict.items()}
        #         elem = QuantumDoubleSchubertAlgebraElement(elem_dict, self)
        #         if cv is not None:
        #             elem = self([1, 2], cv) * elem
        #     else:
        #         return self(x.expand(), cv)
        elif isinstance(x, QuantumDoubleSchubertAlgebraElement):
            if x.is_Add or x.is_Mul:
                return x
            if x._base_var == self._base_var:
                elem = QuantumDoubleSchubertAlgebraElement(x._doubledict)  # , self)
            else:
                return self(x.expand(), cv)
        else:
            x = sympify(x)
            if cv is None or cv == utils.NoneVar:
                cv = utils.NoneVar
                result = py.mult_poly_py({Permutation([]): 1}, x, utils.poly_ring(self._base_var))
            else:
                result = yz.mult_poly_double({Permutation([]): 1}, x, utils.poly_ring(self._base_var), utils.poly_ring(cv))
            elem = QuantumDoubleSchubertAlgebraElement({(k, cv): v for k, v in result.items()})
        return elem


# def _do_schub_mul(a, b):
#     A = QDSx(a)
#     B = QDSx(b)
#     return _from_double_dict(_mul_schub_dicts(A._doubledict, B._doubledict))


# def _do_schub_add(a, b):
#     A = QDSx(a)
#     B = QDSx(b)
#     return _from_double_dict(add_perm_dict(A._doubledict, B._doubledict))


# def get_postprocessor(cls):
#     if cls is Mul:
#         return lambda expr: SchubMul(*expr.args, evaluate=False)  # .doit()
#     if cls is Add:
#         return lambda expr: SchubAdd(*expr.args, evaluate=False)  # .doit()
#     return None

    


# # Basic._constructor_postprocessor_mapping[DoubleSchubertAlgebraElement] = {
# #     "Mul": [get_postprocessor(Mul)],
# #     "Add": [get_postprocessor(Add)],
# # }

# # add.register_handlerclass((Expr, SchubAdd), SchubAdd)
# # mul.register_handlerclass((Expr, SchubMul), SchubMul)


QDSx = QuantumDoubleSchubertAlgebraElement_basis()


def QSx(x):
    return QDSx(x, utils.NoneVar)


QuantumDoubleSchubertPolynomial = QuantumDoubleSchubertAlgebraElement

# is_Add = True
# is_Mul = True
# is_Add
# is_AlgebraicNumber
# is_Atom
# is_Boolean
# is_Derivative
# is_Dummy
# is_Equality
# is_Float
# is_Function
# is_Indexed
# is_Integer
# is_MatAdd
# is_MatMul
# is_Matrix
# is_Mul
# is_Not
# is_Number
# is_NumberSymbol
# is_Order
# is_Piecewise
# is_Point
# is_Poly
# is_Pow
# is_Rational
# is_Relational
# is_Symbol
# is_Vector
# is_Wild
# is_algebraic
# is_algebraic_expr
# is_antihermitian
# is_commutative
# is_comparable
# is_complex
# is_composite
# is_constant
# is_even
# is_extended_negative
# is_extended_nonnegative
# is_extended_nonpositive
# is_extended_nonzero
# is_extended_positive
# is_extended_real
# is_finite
# is_hermitian
# is_hypergeometric
# is_imaginary
# is_infinite
# is_integer
# is_irrational
# is_meromorphic
# is_negative
# is_noninteger
# is_nonnegative
# is_nonpositive
# is_nonzero
# is_number
# is_odd
# is_polar
# is_polynomial
# is_positive
# is_prime
# is_rational
# is_rational_function
# is_real
# is_scalar
# is_symbol
# is_transcendental
# is_zero
# is_polynomial = True
# is_Symbol = True


# class SchubAdd(QuantumDoubleSchubertAlgebraElement, Add):
#     is_Add = True

#     def __new__(cls, *args, evaluate=True, _sympify=True):
#         obj = Add.__new__(cls, *args, evaluate=evaluate, _sympify=_sympify)
#         if evaluate:
#             return obj.doit()
#         return obj

#     def doit(self):
#         ret = self.args[0]
#         for arg in self.args[1:]:
#             if arg.is_Add or arg.is_Mul:
#                 arg = arg.doit()
#             ret = _do_schub_add(ret, arg)
#         return ret

#     # def _sympystr(self, printer):
#     #     return _def_printer._print(f"SchubAdd({self.args}")


# class SchubMul(QuantumDoubleSchubertAlgebraElement, Mul):
#     is_Mul = True

#     def __new__(cls, *args, evaluate=True, _sympify=True):
#         if len(args) == 0:
#             return 1
#         # args, a, b = Mul.flatten(list(args))
#         # if len(args) == 0:
#         #     return 1
#         obj = Mul.__new__(cls, *args, evaluate=evaluate, _sympify=_sympify)

#         if evaluate:
#             return obj.doit()
#         return obj

#     def doit(self):
#         ret = self.args[0]
#         for arg in self.args[1:]:
#             if arg.is_Add or arg.is_Mul:
#                 arg = arg.doit()
#             ret = _do_schub_mul(ret, arg)
#         return ret


# Basic._constructor_postprocessor_mapping[DoubleSchubertAlgebraElement] = {
#     "Mul": [get_postprocessor(Mul)],
#     "Add": [get_postprocessor(Add)],
# }


