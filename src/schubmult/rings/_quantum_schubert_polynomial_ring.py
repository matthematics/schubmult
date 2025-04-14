# to encourage development

from functools import cache

import numpy as np
import sympy
from symengine import S, sympify
from sympy import Basic

import schubmult.rings._schubert_polynomial_ring as spr
import schubmult.rings._utils as utils
import schubmult.schub_lib.quantum as py
import schubmult.schub_lib.quantum_double as yz
from schubmult.perm_lib.perm_lib import Permutation, count_less_than, is_parabolic, longest_element, omega, permtrim, trimcode
from schubmult.poly_lib.poly_lib import elem_sym_poly_q, q_vector, xreplace_genvars
from schubmult.poly_lib.schub_poly import schubpoly_from_elems
from schubmult.poly_lib.variables import GeneratingSet, GeneratingSet_base
from schubmult.schub_lib.schub_lib import check_blocks
from schubmult.utils.logging import get_logger

## EMULATE POLYTOOLS

q_var = GeneratingSet("q")
# _def_printer = StrPrinter({"order": "none"})

logger = get_logger(__name__)


class QuantumDoubleSchubertAlgebraElement(spr.BasisSchubertAlgebraElement):
    def __new__(cls, _dict, basis):
        return spr.BasisSchubertAlgebraElement.__new__(cls, _dict, basis)

    def subs(self, old, new):
        logger.debug("ferefef")
        elb = self.as_classical().subs(old, new).as_quantum()
        logger.debug(f"{elb=}")
        return elb


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


class ParabolicQuantumDoubleSchubertAlgebraElement(spr.BasisSchubertAlgebraElement):
    def __new__(cls, _dict, basis):
        obj = spr.BasisSchubertAlgebraElement.__new__(cls, _dict, basis)
        # obj._index_comp = tuple(index_comp)parabolic_index = []
        # start = 0
        # # 1, 2 | 3 
        # for i in range(len(args.parabolic)):
        #     end = start + int(args.parabolic[i])
        #     parabolic_index += list(range(start+1,end))
        #     # start += int(args.parabolic[i])
        #     start = end
        return obj

    @property
    def index_comp(self):
        return self.basis.index_comp

    # return (sympy.Dict(self._dict), self._basis)     return obj

    # @property
    # def args(self):
    #     return

    # def _hashable_content(self):
    #     return self.args


class PQDSchubPoly(ParabolicQuantumDoubleSchubertAlgebraElement):
    is_Atom = True

    def __new__(cls, k, basis):
        return PQDSchubPoly.__xnew_cached__(cls, k, basis)

    @staticmethod
    def __xnew__(_class, k, basis):
        obj = ParabolicQuantumDoubleSchubertAlgebraElement.__new__(_class, sympy.Dict({(Permutation(k[0]), k[1]): 1}), basis)
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
            return printer.doprint(f"QDS{self.genset.label}({(list(trimcode(self._perm)), list(self.index_comp))}")
        return printer.doprint(f"QDS{self.genset.label}({(list(trimcode(self._perm)), list(self.index_comp))}, {spr._varstr(self._coeff_var)})")


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
        result = S.Zero
        for k, v in elem.coeff_dict.items():
            result += v * self.quantum_as_classical_schubpoly(k[0], k[1])
        return result

    def classical_elem_func(self, coeff_var):
        basis = spr.DoubleSchubertAlgebraElement_basis(self.genset)
        q_var = yz._vars.q_var

        def elem_func(p, k, varl1, varl2):
            if p == 0 and k >= 0:
                return basis([], coeff_var)
            if p < 0 or p > k:
                return basis(0, coeff_var)
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
        return schubpoly_from_elems(k[0], self.genset, utils.poly_ring(k[1]), elem_func=elem_sym_poly_q)  # yz.schubpoly_quantum(k[0], self.genset, utils.poly_ring(k[1]))

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
        # logger.debug(f"{x=} {type(x)=}")
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
                elem = QuantumDoubleSchubertAlgebraElement(x.coeff_dict, self)  # , self)
            else:
                return self(x.expand(), cv, genset)
        elif isinstance(x, spr.DoubleSchubertAlgebraElement):
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
            elem = QuantumDoubleSchubertAlgebraElement({(k, cv): v for k, v in result.items()}, self)
        return elem


QDSx = QuantumDoubleSchubertAlgebraElement_basis(GeneratingSet("x"))


class QuantumSchubertAlgebraElement_basis(QuantumDoubleSchubertAlgebraElement_basis):
    def __new__(cls, genset):
        return QuantumDoubleSchubertAlgebraElement_basis.__new__(cls, genset)

    def _from_single_dict(self, _dict):
        return QuantumDoubleSchubertAlgebraElement({(k, utils.NoneVar): v for k, v in _dict.items()}, self)

    def __call__(self, x):
        genset = self.genset
        # logger.debug(f"{x=} {type(x)=}")
        if not genset:
            genset = self.genset
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        if isinstance(x, list) or isinstance(x, tuple):
            elem = self._from_single_dict({Permutation(x): 1})
        elif isinstance(x, Permutation):
            elem = self._from_single_dict({x: 1})
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
                elem = QuantumDoubleSchubertAlgebraElement(x.coeff_dict, self)  # , self)
            else:
                return self(x.expand())
        elif isinstance(x, spr.DoubleSchubertAlgebraElement):
            if x.genset == self.genset:
                return x.as_quantum()
        else:
            x = sympify(x)
            result = py.mult_poly_q({Permutation([]): 1}, x, genset)
            elem = self._from_single_dict(result)
        return elem


class ParabolicQuantumDoubleSchubertAlgebraElement_basis(Basic):
    def __new__(cls, genset, index_comp):
        obj = Basic.__new__(cls, genset, tuple(index_comp))
        parabolic_index = []
        start = 0
        # 1, 2 | 3 
        for i in range(len(index_comp)):
            end = start + index_comp[i]
            parabolic_index += list(range(start+1,end))
            # start += int(args.parabolic[i])
            start = end
        obj._parabolic_index = parabolic_index
        return obj
        

    def _from_dict(self, _dict):
        return ParabolicQuantumDoubleSchubertAlgebraElement(_dict, self)

    @property
    def genset(self):
        return self.args[0]

    @property
    def index_comp(self):
        return self.args[1]

    def process_coeff_dict(self, coeff_dict):
        max_len = max(len(w) for w in coeff_dict)
        parabolic_index = [*self._parabolic_index]
        parabolic_index += list(range(parabolic_index[-1] + 2, max_len))
        w_P = longest_element(parabolic_index)
        # max_len = len(w_P)
        w_P_prime = Permutation([1, 2])
        coeff_dict_update = {}
        for w_1 in coeff_dict.keys():
            val = coeff_dict[w_1]
            q_dict = yz.factor_out_q_keep_factored(val)
            for q_part in q_dict:
                qv = q_vector(q_part)
                w = w_1
                good = True
                parabolic_index2 = []
                for i in range(len(parabolic_index)):
                    if omega(parabolic_index[i], qv) == 0:
                        parabolic_index2 += [parabolic_index[i]]
                    elif omega(parabolic_index[i], qv) != -1:
                        good = False
                        break
                if not good:
                    continue
                w_P_prime = longest_element(parabolic_index2)
                if not check_blocks(qv, parabolic_index):
                    continue
                w = (w * w_P_prime) * w_P
                if not is_parabolic(w, parabolic_index):
                    continue

                w = permtrim(w)
                if len(w) > max_len:
                    continue
                new_q_part = np.prod(
                    [q_var[index + 1 - count_less_than(parabolic_index, index + 1)] ** qv[index] for index in range(len(qv)) if index + 1 not in parabolic_index],
                )
                try:
                    new_q_part = int(new_q_part)
                except Exception:
                    pass
                q_val_part = q_dict[q_part]
                coeff_dict_update[w] = coeff_dict_update.get(w, 0) + new_q_part * q_val_part
        return coeff_dict_update

    @cache
    def cached_product(self, u, v, va, vb):
        initial_dict = {k: xreplace_genvars(x, utils.poly_ring(va), utils.poly_ring(vb)) for k, x in yz.schubmult_q_double_pair_generic(u, v).items()}
        return {(k, va): v for k, v in self.process_coeff_dict(initial_dict).items()}

    def in_quantum_basis(self, elem):
        return elem

    def in_classical_basis(self, elem):
        result = S.Zero
        for k, v in elem.coeff_dict.items():
            result += v * self.quantum_as_classical_schubpoly(k[0], k[1])
        return result

    def classical_elem_func(self, coeff_var):
        basis = spr.DoubleSchubertAlgebraElement_basis(self.genset)
        q_var = yz._vars.q_var

        def elem_func(p, k, varl1, varl2):
            if p == 0 and k >= 0:
                return basis([], coeff_var)
            if p < 0 or p > k:
                return basis(0, coeff_var)
            return (varl1[k - 1] - varl2[k - p]) * elem_func(p - 1, k - 1, varl1, varl2) + elem_func(p, k - 1, varl1, varl2) + q_var[k - 1] * elem_func(p - 2, k - 2, varl1, varl2)

        return elem_func

    @property
    def single_element_class(self):
        return PQDSchubPoly

    @cache
    def quantum_as_classical_schubpoly(self, perm, coeff_var="y"):
        return schubpoly_from_elems(perm, self.genset, utils.poly_ring(coeff_var), self.classical_elem_func(coeff_var))

    @cache
    def cached_schubpoly(self, k):
        return schubpoly_from_elems(k[0], self.genset, utils.poly_ring(k[1]), elem_func=elem_sym_poly_q)  # yz.schubpoly_quantum(k[0], self.genset, utils.poly_ring(k[1]))

    @cache
    def cached_positive_product(self, u, v, va, vb):
        initial_dict = {k: xreplace_genvars(x, utils.poly_ring(va), utils.poly_ring(vb)) for k, x in yz.schubmult_q_generic_partial_posify(u, v).items()}
        return {(k, va): v for k, v in self.process_coeff_dict(initial_dict).items()}

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
        # logger.debug(f"{x=} {type(x)=}")
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
        elif isinstance(x, ParabolicQuantumDoubleSchubertAlgebraElement):
            return x
        else:
            raise NotImplementedError
        # elif isinstance(x, spr.SchubertPolynomial):
        #     if x._parent._base_var == self._base_var:
        #         elem_dict = {(x, utils.NoneVar): v for k, v in x.coeff_dict.items()}
        #         elem = QuantumDoubleSchubertAlgebraElement(elem_dict, self)
        #         if cv is not None:
        #             elem = self([1, 2], cv) * elem
        #     else:
        #         return self(x.expand(), cv)
        # elif isinstance(x, QuantumDoubleSchubertAlgebraElement):
        #     if x.is_Add or x.is_Mul:
        #         return x
        #     if x.genset == genset:
        #         elem = QuantumDoubleSchubertAlgebraElement(x.coeff_dict, self)  # , self)
        #     else:
        #         return self(x.expand(), cv, genset)
        # elif isinstance(x, spr.DoubleSchubertAlgebraElement):
        #     if x.genset == self.genset:
        #         return self(x.expand(), cv, genset)
        # else:
        #     logger.debug("bagelflap")
        #     x = sympify(x)
        #     if cv is None or cv == utils.NoneVar:
        #         cv = utils.NoneVar
        #         logger.debug(f"{x=} {list(genset)=}")
        #         result = py.mult_poly_q({Permutation([]): 1}, x, genset)
        #         logger.debug(f"{result=}")
        #     else:
        #         result = yz.mult_poly_q_double({Permutation([]): 1}, x, genset, utils.poly_ring(cv))
        #     elem = QuantumDoubleSchubertAlgebraElement({(k, cv): v for k, v in result.items()}, self)
        return elem


QSx = QuantumSchubertAlgebraElement_basis(GeneratingSet("x"))

QuantumDoubleSchubertPolynomial = QuantumDoubleSchubertAlgebraElement


def make_parabolic_quantum_basis(index_comp):
    return ParabolicQuantumDoubleSchubertAlgebraElement_basis(GeneratingSet("x"), index_comp)


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
