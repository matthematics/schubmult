# to encourage development

from functools import cache

import sympy
from symengine import expand, sympify
from sympy import Basic

import schubmult.rings._schubert_polynomial_ring as spr
import schubmult.rings._utils as utils
import schubmult.schub_lib.quantum as py
import schubmult.schub_lib.quantum_double as yz
from schubmult.perm_lib import (
    Permutation,
)
from schubmult.poly_lib.poly_lib import elem_sym_poly_q, xreplace_genvars
from schubmult.poly_lib.schub_poly import schubpoly_from_elems
from schubmult.poly_lib.variables import GeneratingSet, GeneratingSet_base
from schubmult.rings._schubert_polynomial_ring import DoubleSchubertAlgebraElement
from schubmult.utils.logging import get_logger

q_var = GeneratingSet("q")


## EMULATE POLYTOOLS


# _def_printer = StrPrinter({"order": "none"})

logger = get_logger(__name__)


class QuantumDoubleSchubertAlgebraElement(spr.BasisSchubertAlgebraElement):
    def __new__(cls, _dict, basis):
        _dict = {k: v for k, v in _dict.items() if expand(v) != 0}
        return QuantumDoubleSchubertAlgebraElement.__xnew_cached__(cls, sympy.Dict(_dict), basis)

    @staticmethod
    def __xnew__(_class, _dict, basis):
        return spr.BasisSchubertAlgebraElement.__new__(_class, _dict, basis)

    @staticmethod
    @cache
    def __xnew_cached__(_class, _dict, basis):
        return QuantumDoubleSchubertAlgebraElement.__xnew__(_class, _dict, basis)

    def _eval_subs(self, old, new):
        return self.as_classical().subs(old, new).as_quantum()


# # TODO: not a noncommutative symbol, something else
# # Atomic Schubert polynomial
class QDSchubPoly(QuantumDoubleSchubertAlgebraElement):
    is_Atom = True

    def __new__(cls, k, basis):
        return QDSchubPoly.__xnew_cached__(cls, k, basis)

    @staticmethod
    def __xnew__(_class, k, basis):
        obj = QuantumDoubleSchubertAlgebraElement.__new__(_class, sympy.Dict({(Permutation(k[0]), k[1]): 1}), basis)
        obj._perm = k[0]
        obj._coeff_var = k[1]
        # obj._base_var = base_var
        return obj

    @staticmethod
    @cache
    def __xnew_cached__(_class, k, genset):
        return QDSchubPoly.__xnew__(_class, k, genset)

    def _sympystr(self, printer):
        if self._coeff_var == 0 or self._coeff_var == utils.NoneVar:
            return printer.doprint(f"QS{self.genset.label}({list(self._perm)})")
        return printer.doprint(f"QDS{self.genset.label}({list(self._perm)}, {spr._varstr(self._coeff_var)})")


# None is faster to store
# class QuantumDoubleSchubertAlgebraElement_basis(Basic):
#     def __new__(cls, genset):
#         return Basic.__new__(cls, genset)

#     def _from_dict(self, _dict):
#         return QuantumDoubleSchubertAlgebraElement(_dict, self.genset)

#     @property
#     def genset(self):
#         return self.args[0]

#     def __call__(self, x, cv=None, genset=None):
#         logger.debug(f"{x=} {type(x)=}")
#         if not genset:
#             genset = self.genset
#         if not isinstance(genset, GeneratingSet_base):
#             raise TypeError
#         if isinstance(x, list) or isinstance(x, tuple):
#             if cv is None:
#                 cv = "y"
#             elem = QuantumDoubleSchubertAlgebraElement({(Permutation(x), cv): 1}, genset)
#         elif isinstance(x, Permutation):
#             if cv is None:
#                 cv = "y"
#             elem = QuantumDoubleSchubertAlgebraElement({(x, cv): 1}, genset)
#         # elif isinstance(x, spr.SchubertPolynomial):
#         #     if x._parent._base_var == self._base_var:
#         #         elem_dict = {(x, utils.NoneVar): v for k, v in x.coeff_dict.items()}
#         #         elem = QuantumDoubleSchubertAlgebraElement(elem_dict, self)
#         #         if cv is not None:
#         #             elem = self([1, 2], cv) * elem
#         #     else:
#         #         return self(x.expand(), cv)
#         elif isinstance(x, QuantumDoubleSchubertAlgebraElement):
#             if x.is_Add or x.is_Mul:
#                 return x
#             if x.genset == genset:
#                 elem = QuantumDoubleSchubertAlgebraElement(x.coeff_dict, genset)  # , self)
#             else:
#                 return self(x.expand(), cv, genset)
#         elif isinstance(x, DoubleSchubertAlgebraElement):
#             if x.genset == self.genset:
#                 return self(x.expand(), cv, genset)
#         else:
#             logger.debug("bagelflap")
#             x = sympify(x)
#             if cv is None or cv == utils.NoneVar:
#                 cv = utils.NoneVar
#                 logger.debug(f"{x=} {list(genset)=}")
#                 result = py.mult_poly_q({Permutation([]): 1}, x, genset)
#                 logger.debug(f"{result=}")
#             else:
#                 result = yz.mult_poly_q_double({Permutation([]): 1}, x, genset, utils.poly_ring(cv))
#             elem = QuantumDoubleSchubertAlgebraElement({(k, cv): v for k, v in result.items()}, genset)
#         return elem


# def _do_schub_mul(a, b):
#     A = QDSx(a)
#     B = QDSx(b)
#     return self._from_dict(_mul_q_schub_dicts(A.coeff_dict, B.coeff_dict))


# def _do_schub_add(a, b):
#     A = QDSx(a)
#     B = QDSx(b)
#     return self._from_dict(add_perm_dict(A.coeff_dict, B.coeff_dict))


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


class QuantumDoubleSchubertAlgebraElement_basis(Basic):
    def __new__(cls, genset):
        return Basic.__new__(cls, genset)

    def _from_dict(self, _dict):
        return QuantumDoubleSchubertAlgebraElement(_dict, self)

    @property
    def genset(self):
        return self.args[0]

    @cache
    def cached_product(self, u, v, va, vb):
        return {(k, va): xreplace_genvars(x, utils.poly_ring(va), utils.poly_ring(vb)) for k, x in yz.schubmult_q_double_pair_generic(u, v).items()}

    def in_quantum_basis(self, elem):
        return elem

    def in_classical_basis(self, elem):
        result = 0
        for k, v in elem.coeff_dict.items():
            result += v * self.quantum_as_classical_schubpoly(k[0], k[1])
        return result

    def classical_elem_func(self, coeff_var):
        basis = spr.DoubleSchubertAlgebraElement_basis(self.genset)
        def elem_func(p, k, varl1, varl2):
            if p == 0 and k >= 0:
                return basis([], coeff_var)
            if p < 0 or p > k:
                return 0
            return (varl1[k - 1] - varl2[k - p]) * elem_func(p - 1, k - 1, varl1, varl2) + elem_func(p, k - 1, varl1, varl2) + q_var[k - 1] * elem_func(p - 2, k - 2, varl1, varl2)

        return elem_func


    @property
    def single_element_class(self):
        return QDSchubPoly

    @cache
    def quantum_as_classical_schubpoly(self, perm, coeff_var="y"):
        return schubpoly_from_elems(perm, self.genset, utils.poly_ring(coeff_var), self.classical_elem_func(coeff_var))

    @cache
    def cached_schubpoly(self, k):
        return schubpoly_from_elems(k[0],self.genset,utils.poly_ring(k[1]),elem_func=elem_sym_poly_q)# yz.schubpoly_quantum(k[0], self.genset, utils.poly_ring(k[1]))

    @cache
    def cached_positive_product(self, u, v, va, vb):
        return {(k, va): xreplace_genvars(x, utils.poly_ring(va), utils.poly_ring(vb)) for k, x in yz.schubmult_q_generic_partial_posify(u, v).items()}

    @property
    def double_mul(self):
        return yz.schubmult_q_double

    @property
    def single_mul(self):
        return py.schubmult_q_fast

    @property
    def mult_poly_single(self):
        return py.mult_poly_q

    @property
    def mult_poly_double(self):
        return yz.mult_poly_q_double

    def __call__(self, x, cv=None):
        genset = self.genset
        logger.debug(f"{x=} {type(x)=}")
        if not genset:
            genset = self.genset
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        if isinstance(x, list) or isinstance(x, tuple):
            if cv is None:
                cv = "y"
            elem = self._from_dict({(Permutation(x), cv): 1})
        elif isinstance(x, Permutation):
            if cv is None:
                cv = "y"
            elem = self._from_dict({(x, cv): 1})
        # elif isinstance(x, spr.SchubertPolynomial):
        #     if x._parent._base_var == self._base_var:
        #         elem_dict = {(x, utils.NoneVar): v for k, v in x.coeff_dict.items()}
        #         elem = QuantumDoubleSchubertAlgebraElement(elem_dict, self)
        #         if cv is not None:
        #             elem = self([1, 2], cv) * elem
        #     else:
        #         return self(x.expand(), cv)
        elif isinstance(x, QuantumDoubleSchubertAlgebraElement):
            if x.is_Add or x.is_Mul:
                return x
            if x.genset == genset:
                elem = QuantumDoubleSchubertAlgebraElement(x.coeff_dict, genset)  # , self)
            else:
                return self(x.expand(), cv, genset)
        elif isinstance(x, DoubleSchubertAlgebraElement):
            if x.genset == self.genset:
                return self(x.expand(), cv, genset)
        else:
            logger.debug("bagelflap")
            x = sympify(x)
            if cv is None or cv == utils.NoneVar:
                cv = utils.NoneVar
                logger.debug(f"{x=} {list(genset)=}")
                result = py.mult_poly_q({Permutation([]): 1}, x, genset)
                logger.debug(f"{result=}")
            else:
                result = yz.mult_poly_q_double({Permutation([]): 1}, x, genset, utils.poly_ring(cv))
            elem = QuantumDoubleSchubertAlgebraElement({(k, cv): v for k, v in result.items()}, genset)
        return elem


QDSx = QuantumDoubleSchubertAlgebraElement_basis(GeneratingSet("x"))


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
