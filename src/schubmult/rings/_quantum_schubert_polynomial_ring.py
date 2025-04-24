# to encourage development

from bisect import bisect_left
from functools import cache

import sympy
from symengine import S, Symbol, sympify
from sympy import Basic

import schubmult.rings._schubert_polynomial_ring as spr
import schubmult.rings._utils as utils
import schubmult.schub_lib.quantum as py
import schubmult.schub_lib.quantum_double as yz
from schubmult.perm_lib import Permutation, longest_element
from schubmult.poly_lib.poly_lib import complete_sym_poly, elem_sym_poly, elem_sym_poly_q, xreplace_genvars
from schubmult.poly_lib.schub_poly import schubpoly_from_elems
from schubmult.poly_lib.variables import GeneratingSet, GeneratingSet_base
from schubmult.utils.logging import get_logger
from schubmult.utils.perm_utils import is_parabolic

## EMULATE POLYTOOLS

q_var = GeneratingSet("q")
# _def_printer = StrPrinter({"order": "none"})

logger = get_logger(__name__)


class QuantumDoubleSchubertAlgebraElement(spr.BasisSchubertAlgebraElement):
    def __new__(cls, _dict, basis):
        return spr.BasisSchubertAlgebraElement.__new__(cls, _dict, basis)

    # def expand(self, *args, **kwargs):
    #     # print("Frofulating bagel")
    #     return super().expand()

    def subs(self, old, new):
        return self.as_classical().subs(old, new).as_quantum()


# # TODO: not a noncommutative symbol, something else
# # Atomic Schubert polynomial
_pretty_schub_char = "𝕼𝔖"  # noqa: RUF001


class QDSchubPoly(QuantumDoubleSchubertAlgebraElement):
    is_Atom = True

    def __new__(cls, k, basis):
        return QDSchubPoly.__xnew_cached__(cls, k, basis)

    @staticmethod
    def __xnew__(_class, k, basis):
        obj = QuantumDoubleSchubertAlgebraElement.__new__(_class, sympy.Dict({Permutation(k): 1}), basis)
        obj._perm = k
        obj._key = k
        # obj._base_var = base_var
        return obj

    @staticmethod
    @cache
    def __xnew_cached__(_class, k, genset):
        return QDSchubPoly.__xnew__(_class, k, genset)

    def _sympystr(self, printer):
        if self.basis.coeff_genset.label is None:
            return printer.doprint(f"QS{self.genset.label}({printer.doprint(self._perm)})")
        return printer.doprint(f"QDS{self.genset.label}({printer.doprint(self._perm)}, {self.basis.coeff_genset.label})")

    def _pretty(self, printer):
        if self._key == Permutation([]):
            return printer._print(1)
        subscript = printer.doprint(int("".join([str(i) for i in self._key])))
        if self.basis.coeff_genset.label is None:
            return printer._print_Function(sympy.Function(f"{_pretty_schub_char}_{subscript}")(sympy.Symbol(self.genset.label)))
        return printer._print_Function(sympy.Function(f"{_pretty_schub_char}_{subscript}")(sympy.Symbol(f"{self.genset.label}; {self.basis.coeff_genset.label}")))

    def _latex(self, printer):
        if self._key == Permutation([]):
            return printer._print(1)
        # subscript = printer._print(int("".join([str(i) for i in self._key])))
        subscript = sympy.sstr(self._key)
        if self.basis.coeff_genset.label is None:
            return printer._print_Function(sympy.Function("\\widetilde{\\mathfrak{S}}" + f"_{'{' + subscript + '}'}")(sympy.Symbol(self.genset.label)))
        return printer._print_Function(sympy.Function("\\widetilde{\\mathfrak{S}}" + f"_{'{' + subscript + '}'}")(sympy.Symbol(f"{self.genset.label}; {self.basis.coeff_genset.label}")))


class ParabolicQuantumDoubleSchubertAlgebraElement(spr.BasisSchubertAlgebraElement):
    def __new__(cls, _dict, basis):
        return spr.BasisSchubertAlgebraElement.__new__(cls, _dict, basis)

    @property
    def index_comp(self):
        return self.basis.index_comp

    def kill_ideal(self):
        length = sum(self.index_comp)
        new_dict = {}
        for k, v in self.coeff_dict.items():
            if len(k) <= length:
                new_dict[k] = v
        return self.basis._from_dict(new_dict)


class PQDSchubPoly(ParabolicQuantumDoubleSchubertAlgebraElement):
    is_Atom = True

    def __new__(cls, k, basis):
        return PQDSchubPoly.__xnew_cached__(cls, k, basis)

    @staticmethod
    def __xnew__(_class, k, basis):
        obj = ParabolicQuantumDoubleSchubertAlgebraElement.__new__(_class, sympy.Dict({Permutation(k): 1}), basis)
        obj._perm = k
        obj._key = k
        # obj._perm._print_as_code = True
        # obj._base_var = base_var
        return obj

    @staticmethod
    @cache
    def __xnew_cached__(_class, k, genset):
        return PQDSchubPoly.__xnew__(_class, k, genset)

    def _sympystr(self, printer):
        if self._coeff_var == 0 or self._coeff_var == utils.NoneVar:
            return printer.doprint(f"QPS{self.genset.label}{(tuple(self.index_comp))}({printer.doprint(self._perm)})")
        return printer.doprint(f"QPDS{self.genset.label}{tuple(self.index_comp)}({printer.doprint(self._perm)}, {spr._varstr(self._coeff_var)})")

    def _pretty(self, printer):
        if self._key == Permutation([]):
            return printer._print(1)
        subscript = printer._print(int("".join([str(i) for i in self._key])))
        # subscript = sympy.sstr(self._key)
        if self._key[1] == 0 or self._key[1] == utils.NoneVar:
            return printer._print_Function(sympy.Function(f"{_pretty_schub_char}_{subscript}")(sympy.Symbol(f"{self.genset.label} | {self.basis.index_comp}")))
        return printer._print_Function(sympy.Function(f"{_pretty_schub_char}_{subscript}")(sympy.Symbol(f"{self.genset.label}; {self._key[1]} | {self.basis.index_comp}")))

    def _latex(self, printer):
        if self._key == Permutation([]):
            return printer._print(1)
        # subscript = printer._print(int("".join([str(i) for i in self._key])))
        subscript = printer._print(self._key)
        supscript = printer._print(sympy.Tuple(*self.index_comp))
        if self._key[1] == 0 or self._key[1] == utils.NoneVar:
            return printer._print_Function(sympy.Function("\\widetilde{\\mathfrak{S}}" + f"^{'{'}{supscript}{'}'}_{'{' + subscript + '}'}")(sympy.Symbol(self.genset.label)))
        return printer._print_Function(sympy.Function("\\widetilde{\\mathfrak{S}}" + f"^{'{'}{supscript}{'}'}_{'{' + subscript + '}'}")(sympy.Symbol(f"{self.genset.label}; {self._key[1]}")))


class QuantumDoubleSchubertAlgebraElement_basis(Basic):
    def __new__(cls, genset, coeff_genset):
        # print(f"{genset=} {coeff_genset=}")
        return QuantumDoubleSchubertAlgebraElement_basis.__xnew_cached__(cls, genset, coeff_genset)

    @staticmethod
    @cache
    def __xnew_cached__(_class, genset, coeff_genset):
        return QuantumDoubleSchubertAlgebraElement_basis.__xnew__(_class, genset, coeff_genset)

    @staticmethod
    def __xnew__(_class, genset, coeff_genset):
        return Basic.__new__(_class, genset, coeff_genset)

    def _coerce_mul(self, other):
        if isinstance(other, spr.BasisSchubertAlgebraElement):
            if type(other.basis) is type(self):
                if self.genset == other.basis.genset:
                    return other
            if type(other.basis) is QuantumSchubertAlgebraElement_basis:
                if self.genset == other.basis.genset:
                    newbasis = QuantumDoubleSchubertAlgebraElement_basis(self.genset, utils.poly_ring(0))
                    return newbasis._from_dict(other.coeff_dict)
        return None

    def _coerce_add(self, other):
        if isinstance(other, spr.BasisSchubertAlgebraElement):
            if type(other.basis) is type(self):
                if self.genset == other.basis.genset and self.coeff_genset == other.basis.coeff_genset:
                    return other
        return None

    @property
    def coeff_genset(self):
        return self.args[1]

    @property
    def symbol_elem_func(self):
        def elem_func(p, k, varl1, varl2):  # noqa: ARG001
            if p == 0 and k >= 0:
                return 1
            if p < 0 or p > k:
                return 0
            return sympy.Add(*[(Symbol(f"e_{p - i}_{k}") if p - i > 0 else 1) * complete_sym_poly(i, k + 1 - p, [-v for v in varl2]) for i in range(p + 1)])

        return elem_func

    def elem_sym_subs(self, kk):
        elems = []
        for k in range(1, kk + 1):
            for p in range(1, k + 1):
                elems += [(sympy.Symbol(f"e_{p}_{k}"), elem_sym_poly_q(p, k, self.genset[1:], utils.poly_ring(0)))]
        return dict(elems)

    def _from_dict(self, _dict):
        return QuantumDoubleSchubertAlgebraElement(_dict, self)

    @property
    def genset(self):
        return self.args[0]

    @cache
    def cached_product(self, u, v, basis2):
        return {k: xreplace_genvars(x, self.coeff_genset, basis2.coeff_genset) for k, x in yz.schubmult_q_double_pair_generic(u, v).items()}

    def in_quantum_basis(self, elem):
        return elem

    def in_classical_basis(self, elem):
        result = S.Zero
        # print(f"{elem=}")
        for k, v in elem.coeff_dict.items():
            result += v * self.quantum_as_classical_schubpoly(k)
            # print(f"{v=} {result=}")
        return result

    @property
    def classical_elem_func(self):
        basis = spr.DoubleSchubertAlgebraElement_basis(self.genset, self.coeff_genset)
        # print(f"{basis=} {self=} {basis.genset=} {self.coeff_genset=} {basis.coeff_genset=}")
        q_var = yz._vars.q_var

        def elem_func(p, k, varl1, varl2):
            if p == 0 and k >= 0:
                return basis([])
            if p < 0 or p > k:
                return basis(0)
            return (varl1[k - 1] - varl2[k - p]) * elem_func(p - 1, k - 1, varl1, varl2) + elem_func(p, k - 1, varl1, varl2) + q_var[k - 1] * elem_func(p - 2, k - 2, varl1, varl2)

        return elem_func

    @property
    def single_element_class(self):
        return QDSchubPoly

    @cache
    def quantum_as_classical_schubpoly(self, perm):
        # print(f"{self=} {self.genset=} {self.coeff_genset=} {perm=}")
        return schubpoly_from_elems(perm, self.genset, self.coeff_genset, self.classical_elem_func)

    @cache
    def cached_schubpoly(self, k):
        # print(f"pralfasank {self=}")
        return schubpoly_from_elems(k, self.genset, self.coeff_genset, elem_func=elem_sym_poly_q)

    @cache
    def cached_positive_product(self, u, v, basis2):
        return {k: xreplace_genvars(x, self.coeff_genset, basis2.coeff_genset) for k, x in yz.schubmult_q_generic_partial_posify(u, v).items()}

    @property
    def double_mul(self):
        return yz.schubmult_q_double_fast

    @property
    def single_mul(self):
        return py.schubmult_q_fast

    @property
    def mult_poly_single(self):
        return py.mult_poly_q

    @property
    def mult_poly_double(self):
        return yz.mult_poly_q_double

    def __call__(self, x):
        genset = self.genset
        # logger.debug(f"{x=} {type(x)=}")
        if not genset:
            genset = self.genset
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        if isinstance(x, list) or isinstance(x, tuple):
            elem = self._from_dict({Permutation(x): 1})
        elif isinstance(x, Permutation):
            elem = self._from_dict({x: 1})
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
        else:
            # logger.debug("bagelflap")
            x = sympify(x)
            result = yz.mult_poly_q_double({Permutation([]): 1}, x, self.genset, self.coeff_genset)
            elem = QuantumDoubleSchubertAlgebraElement(result, self)
        return elem


def QDSx(x, genset=GeneratingSet("y")):
    if isinstance(genset, str):
        genset = GeneratingSet(genset)
    # print(f"{genset=}")
    return QuantumDoubleSchubertAlgebraElement_basis(GeneratingSet("x"), genset)(x)


t = GeneratingSet("t")
a = GeneratingSet("a")

spunky_basis = spr.SchubertAlgebraElement_basis(t)


class QuantumSchubertAlgebraElement_basis(QuantumDoubleSchubertAlgebraElement_basis):
    def __new__(cls, genset):
        return QuantumDoubleSchubertAlgebraElement_basis.__new__(cls, genset, utils.poly_ring(utils.NoneVar))

    def _coerce_mul(self, other):
        """Coerce a basis schubert algebra element so it can be multiplied

        Args:
            other (_type_): _description_

        Returns:
            _type_: _description_
        """
        if type(other.basis) is type(self):
            if self.genset == other.basis.genset:
                return other
        if type(other.basis) is QuantumDoubleSchubertAlgebraElement_basis:
            if self.genset == other.basis.genset:
                return other
        return None

    @property
    def coeff_genset(self):
        return utils.poly_ring(utils.NoneVar)

    @cache
    def cached_product(self, u, v, basis2):
        if self == basis2:
            return py.schubmult_q_fast({u: S.One}, v)
        return {k: xreplace_genvars(x, self.coeff_genset, basis2.coeff_genset) for k, x in yz.schubmult_q_double_pair_generic(u, v).items()}

    @cache
    def cached_positive_product(self, u, v, basis2):
        return self.cached_product(u, v, basis2)

    def __call__(self, x):
        genset = self.genset
        # logger.debug(f"{x=} {type(x)=}")
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        if isinstance(x, list) or isinstance(x, tuple):
            elem = self._from_dict({Permutation(x): 1})
        elif isinstance(x, Permutation):
            elem = self._from_dict({x: 1})
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
        elif isinstance(x, ParabolicQuantumDoubleSchubertAlgebraElement):
            return x.as_quantum()
        else:
            x = sympify(x)
            result = py.mult_poly_q({Permutation([]): 1}, x, genset)
            elem = self._from_dict(result)
        return elem

    def in_SEM_basis(self):
        result = S.Zero
        for k, v in self.coeff_dict.items():
            if len(k) > len(self._longest):
                parabolic_index = []
                start = 0
                # 1, 2 | 3
                index_comp = [*self._n, len(k) + 1 - self._N[-1]]
                for i in range(len(index_comp)):
                    end = start + index_comp[i]
                    parabolic_index += list(range(start + 1, end))
                    # start += int(args.parabolic[i])
                    start = end
                otherlong = Permutation(list(range(parabolic_index[-1] + 1, 0, -1)))
                longpar = Permutation(longest_element(parabolic_index))
                # print(f"{longpar=} {parabolic_index=}")
                longest = otherlong * longpar
                # print(f"new longest = {longest=}")
            else:
                longest = self._longest
            result += v * schubpoly_from_elems(k, self.genset, self.coeff_genset, elem_func=self.basis.symbol_elem_func, mumu=~longest)
        return result


class ParabolicQuantumDoubleSchubertAlgebraElement_basis(Basic):
    def __new__(cls, genset, coeff_genset, index_comp):
        return ParabolicQuantumDoubleSchubertAlgebraElement_basis.__xnew_cached__(cls, genset, coeff_genset, sympy.Tuple(*index_comp))

    @staticmethod
    @cache
    def __xnew_cached__(_class, genset, coeff_genset, index_comp):
        return ParabolicQuantumDoubleSchubertAlgebraElement_basis.__xnew__(_class, genset, coeff_genset, index_comp)

    @staticmethod
    def __xnew__(_class, genset, coeff_genset, index_comp):
        obj = Basic.__new__(_class, genset, coeff_genset, index_comp)
        obj._quantum_basis = QuantumDoubleSchubertAlgebraElement_basis(genset, coeff_genset)
        obj._classical_basis = spr.DoubleSchubertAlgebraElement_basis(genset, coeff_genset)
        obj._n = list(index_comp)
        obj._N = [sum(obj._n[:i]) for i in range(len(obj._n) + 1)]
        parabolic_index = []
        start = 0
        for i in range(len(index_comp)):
            end = start + index_comp[i]
            parabolic_index += list(range(start + 1, end))
            # start += int(args.parabolic[i])
            start = end
        obj._parabolic_index = parabolic_index
        obj._otherlong = Permutation(list(range(obj._N[-1], 0, -1)))
        obj._longest = obj._otherlong * longest_element(parabolic_index)
        return obj

    def _coerce_mul(self, other):
        if isinstance(other, spr.BasisSchubertAlgebraElement):
            if type(other.basis) is type(self):
                if self.genset == other.basis.genset:
                    return other
            if type(other.basis) is spr.ParabolicQuantumSchubertAlgebraElement_basis:
                if self.genset == other.basis.genset and self.index_comp == other.basis.index_comp:
                    newbasis = ParabolicQuantumDoubleSchubertAlgebraElement_basis(self.genset, utils.poly_ring(0), self.index_comp)
                    return newbasis._from_dict(other.coeff_dict)
        return None

    def _coerce_add(self, other):
        if isinstance(other, spr.BasisSchubertAlgebraElement):
            if type(other.basis) is type(self):
                if self.genset == other.basis.genset and self.coeff_genset == other.basis.coeff_genset and self.index_comp == other.basis.index_comp:
                    return other
        return None

    @property
    def symbol_elem_func(self):
        def elem_func(p, k, varl1, varl2):  # noqa: ARG001
            if p == 0 and k >= 0:
                return 1
            if p < 0 or p > k:
                return 0
            return sympy.Add(*[(Symbol(f"e_{p - i}_{k}") if p - i > 0 else 1) * complete_sym_poly(i, k + 1 - p, [-v for v in varl2]) for i in range(p + 1)])

        return elem_func

    def elem_sym_subs(self, kk):
        elems = []
        elem_func = self.elem_sym
        for k in range(1, kk + 1):
            for p in range(1, k + 1):
                elems += [(sympy.Symbol(f"e_{p}_{k}"), elem_func(p, k, self.genset[1:], utils.poly_ring(0)))]
        return dict(elems)

    @property
    def parabolic_index(self):
        return self._parabolic_index

    @property
    def quantum_basis(self):
        return self._quantum_basis

    @property
    def classical_basis(self):
        return self._classical_basis

    # def elem_sym_poly(self, p, k, varl1, varl2, xstart=0, ystart=0):
    #     # print(f"{p=} {k=} {xstart=} {ystart=} {len(varl1)=} {len(varl2)=}")
    #     if p > k:
    #         return zero
    #     if p == 0:
    #         return one
    #     if p == 1:
    #         res = varl1[xstart] - varl2[ystart]
    #         for i in range(1, k):
    #             res += varl1[xstart + i] - varl2[ystart + i]
    #         return res
    #     if p == k:
    #         res = (varl1[xstart] - varl2[ystart]) * (varl1[xstart + 1] - varl2[ystart])
    #         for i in range(2, k):
    #             res *= varl1[i + xstart] - varl2[ystart]
    #         return res
    #     mid = k // 2
    #     xsm = xstart + mid
    #     ysm = ystart + mid
    #     kmm = k - mid
    #     res = elem_sym_poly(p, mid, varl1, varl2, xstart, ystart) + elem_sym_poly(
    #         p,
    #         kmm,
    #         varl1,
    #         varl2,
    #         xsm,
    #         ysm,
    #     )
    #     for p2 in range(max(1, p - kmm), min(p, mid + 1)):
    #         res += elem_sym_poly(p2, mid, varl1, varl2, xstart, ystart) * elem_sym_poly(
    #             p - p2,
    #             kmm,
    #             varl1,
    #             varl2,
    #             xsm,
    #             ysm - p2,
    #         )
    #     return res

    def elem_sym(self, p, k, varl1, varl2):
        # print(f"{p=} {k=} {self._N=}")
        if p < 0 or p > k:
            return 0
        if p == 0 and k >= 0:
            return 1
        if k <= self._N[1]:
            return elem_sym_poly(p, k, varl1, varl2)
        ret = 0
        j = bisect_left(self._N, k)
        if j < len(self._N) and k == self._N[j]:
            ret = (-((-1) ** (self._n[j - 1]))) * q_var[j - 1] * self.elem_sym(p - self._N[j] + self._N[j - 2], self._N[j - 2], varl1, varl2)
        ret += self.elem_sym(p, k - 1, varl1, varl2) + (varl1[k - 1] - varl2[k - p]) * self.elem_sym(p - 1, k - 1, varl1, varl2)
        return ret

    # def classical_elem(self, k):
    #     if k <= self._N[1]:
    #         return self(uncode([1 for i in range(k)]))

    def _from_dict(self, _dict):
        return ParabolicQuantumDoubleSchubertAlgebraElement(_dict, self)

    @property
    def genset(self):
        return self.args[0]

    @property
    def coeff_genset(self):
        return self.args[1]

    @property
    def index_comp(self):
        return self.args[2]

    def process_coeff_dict(self, coeff_dict):
        max_len = max(len(w) for w in coeff_dict)
        parabolic_index = [*self._parabolic_index]
        # print(f"bagels = {parabolic_index=} {type(self)=}")
        # parabolic_index += list(range(parabolic_index[-1] + 2, max_len + 1))
        if max_len > len(self._longest):
            parabolic_index = []
            start = 0
            # 1, 2 | 3
            index_comp = [*self._n, max_len + 1 - self._N[-1]]
            for i in range(len(index_comp)):
                end = start + index_comp[i]
                parabolic_index += list(range(start + 1, end))
                # start += int(args.parabolic[i])
                start = end
        return yz.apply_peterson_woodward(coeff_dict, parabolic_index)

    @cache
    def cached_product(self, u, v, basis2):
        initial_dict = {k: xreplace_genvars(x, self.coeff_genset, basis2.coeff_genset) for k, x in yz.schubmult_q_double_pair_generic(u, v).items()}
        return self.process_coeff_dict(initial_dict)

    def in_quantum_basis(self, elem):
        result = S.Zero
        for k, v in elem.coeff_dict.items():
            result += v * schubpoly_from_elems(k, self.genset, self.coeff_genset, self.quantum_elem_func)
            # print(f"{result=}")
        return result

    def in_classical_basis(self, elem):
        result = S.Zero
        for k, v in elem.coeff_dict.items():
            result += v * self.quantum_as_classical_schubpoly(k)
            # print(f"{result=}")
        return result

    @cache
    def classical_in_basis(self, k):
        from symengine import expand

        a = self.classical_basis(k)
        b = self(k)
        if expand(a.as_polynomial() - b.as_polynomial()) == S.Zero:
            return b

        cd = dict(b.as_classical().coeff_dict)
        for k2, v in cd.items():
            if k != k2:
                b -= v * self.classical_in_basis(k2)
        return b

    @property
    def classical_elem_func(self):
        basis = spr.DoubleSchubertAlgebraElement_basis(self.genset, self.coeff_genset)
        q_var = yz._vars.q_var

        def elem_func(p, k, varl1, varl2):
            # print(f"{p=} {k=} {varl1=} {varl2=}")
            if p == 0 and k >= 0:
                return basis([])
            if p < 0 or p > k:
                return basis(0)
            if k <= self._N[1]:
                return basis(elem_sym_poly(p, k, varl1, varl2))
            ret = basis(0)
            j = bisect_left(self._N, k)
            if j < len(self._N) and k == self._N[j]:
                ret = (-((-1) ** (self._n[j - 1]))) * q_var[j - 1] * elem_func(p - self._N[j] + self._N[j - 2], self._N[j - 2], varl1, varl2)
            ret += elem_func(p, k - 1, varl1, varl2) + (varl1[k - 1] - varl2[k - p]) * elem_func(p - 1, k - 1, varl1, varl2)
            # print(f"{ret=}")
            return ret

        return elem_func

    @property
    def quantum_elem_func(self):
        basis = QuantumDoubleSchubertAlgebraElement_basis(self.genset, self.coeff_genset)
        q_var = yz._vars.q_var

        def elem_func(p, k, varl1, varl2):
            # print(f"{p=} {k=} {varl1=} {varl2=}")
            if p == 0 and k >= 0:
                return basis([])
            if p < 0 or p > k:
                return basis(0)
            if k <= self._N[1]:
                return basis(elem_sym_poly(p, k, varl1, varl2))
            ret = basis(0)
            j = bisect_left(self._N, k)
            if j < len(self._N) and k == self._N[j]:
                ret = (-((-1) ** (self._n[j - 1]))) * q_var[j - 1] * elem_func(p - self._N[j] + self._N[j - 2], self._N[j - 2], varl1, varl2)
            ret += elem_func(p, k - 1, varl1, varl2) + (varl1[k - 1] - varl2[k - p]) * elem_func(p - 1, k - 1, varl1, varl2)
            # print(f"{ret=}")
            return ret

        return elem_func

    @property
    def single_element_class(self):
        return PQDSchubPoly

    @cache
    def quantum_as_classical_schubpoly(self, perm):
        k = perm
        # print(f"{k=}")
        if len(k) > len(self._longest):
            parabolic_index = []
            start = 0
            # 1, 2 | 3
            index_comp = [*self._n, len(k) + 1 - self._N[-1]]
            for i in range(len(index_comp)):
                end = start + index_comp[i]
                parabolic_index += list(range(start + 1, end))
                # start += int(args.parabolic[i])
                start = end
            otherlong = Permutation(list(range(parabolic_index[-1] + 1, 0, -1)))
            longpar = Permutation(longest_element(parabolic_index))
            # print(f"{longpar=} {parabolic_index=}")
            longest = otherlong * longpar
            # print(f"new longest = {longest=}")
        else:
            longest = self._longest
        return schubpoly_from_elems(perm, self.genset, self.coeff_genset, elem_func=self.classical_elem_func, mumu=~longest)

    @cache
    def cached_schubpoly(self, k):
        if len(k) > len(self._longest):
            parabolic_index = []
            start = 0
            # 1, 2 | 3
            index_comp = [*self._n, len(k) + 1 - self._N[-1]]
            for i in range(len(index_comp)):
                end = start + index_comp[i]
                parabolic_index += list(range(start + 1, end))
                # start += int(args.parabolic[i])
                start = end
            otherlong = Permutation(list(range(parabolic_index[-1] + 1, 0, -1)))
            longpar = Permutation(longest_element(parabolic_index))
            # print(f"{longpar=} {parabolic_index=}")
            longest = otherlong * longpar
            # print(f"new longest = {longest=}")
        else:
            longest = self._longest
        return schubpoly_from_elems(k, self.genset, self.coeff_genset, elem_func=self.elem_sym, mumu=~longest)

    @cache
    def cached_positive_product(self, u, v, basis2):
        initial_dict = {
            k: xreplace_genvars(x, self.coeff_genset, basis2.coeff_genset if basis2.coeff_genset else utils.poly_ring(utils.NoneVar)) for k, x in yz.schubmult_q_generic_partial_posify(u, v).items()
        }
        return self.process_coeff_dict(initial_dict)

    @property
    def double_mul(self):
        from schubmult.schub_lib.quantum_double import _vars

        def do_double_mul(perm_dict, v, var2=_vars.var2, var3=_vars.var3, q_var=_vars.q_var):
            coeff_dict = yz.schubmult_q_double_fast(perm_dict, v, var2, var3, q_var)
            return self.process_coeff_dict(coeff_dict)

        return do_double_mul

    @property
    def single_mul(self):
        from schubmult.schub_lib.quantum_double import _vars

        def do_single_mul(perm_dict, v, q_var=_vars.q_var):
            coeff_dict = py.schubmult_q_fast(perm_dict, v, q_var)
            return self.process_coeff_dict(coeff_dict)

        return do_single_mul

    @property
    def mult_poly_single(self):
        return py.mult_poly_q

    @property
    def mult_poly_double(self):
        return yz.mult_poly_q_double

    def __call__(self, x):
        genset = self.genset
        # logger.debug(f"{x=} {type(x)=}")
        if not genset:
            genset = self.genset
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        if isinstance(x, list) or isinstance(x, tuple):
            perm = Permutation(x)
            if not is_parabolic(perm, self.parabolic_index):
                raise ValueError(f"Permutation must be parabolic: {perm} is not")
            elem = self._from_dict({perm: 1})
        elif isinstance(x, Permutation):
            if not is_parabolic(x, self.parabolic_index):
                raise ValueError(f"Permutation must be parabolic: {x} is not")
            elem = self._from_dict({x: 1})
        elif isinstance(x, ParabolicQuantumDoubleSchubertAlgebraElement):
            return x
        else:
            dct = self.classical_basis(x)
            elem = 0
            if not isinstance(dct, spr.BasisSchubertAlgebraElement):
                return dct
            try:
                for k, v in dct.coeff_dict.items():
                    if elem == 0:
                        elem = v * self.classical_in_basis(k)
                    else:
                        elem += v * self.classical_in_basis(k)
            except ValueError:
                raise ValueError(f"Could not convert {x=} to quantum parabolic")
        return elem


QSx = QuantumSchubertAlgebraElement_basis(GeneratingSet("x"))

QuantumDoubleSchubertPolynomial = QuantumDoubleSchubertAlgebraElement


def make_parabolic_quantum_basis(index_comp, coeff_genset):
    return ParabolicQuantumDoubleSchubertAlgebraElement_basis(GeneratingSet("x"), coeff_genset, index_comp)


def QPDSx_index(*args):
    def this_QPDSx(x, coeff_genset=GeneratingSet("y")):
        return make_parabolic_quantum_basis(args, utils.poly_ring(coeff_genset) if isinstance(coeff_genset, str) else coeff_genset)(x)

    return this_QPDSx


@cache
def QPDSx(*args):
    return QPDSx_index(*args)


class ParabolicQuantumSchubertAlgebraElement_basis(ParabolicQuantumDoubleSchubertAlgebraElement_basis):
    def __new__(cls, genset, index_comp):
        return ParabolicQuantumDoubleSchubertAlgebraElement_basis.__new__(cls, genset, utils.poly_ring(0), index_comp)

    def __hash__(self):
        return hash((*self.args, utils.NoneVar))

    def _coerce_mul(self, other):
        if isinstance(other, spr.BasisSchubertAlgebraElement):
            if type(other.basis) is type(self):
                if self.genset == other.basis.genset:
                    return other
            if type(other.basis) is spr.ParabolicQuantumDoubleSchubertAlgebraElement_basis:
                if self.genset == other.basis.genset and self.index_comp == other.basis.index_comp:
                    newbasis = ParabolicQuantumDoubleSchubertAlgebraElement_basis(self.genset, utils.poly_ring(0), self.index_comp)
                    return newbasis._from_dict(other.coeff_dict)
        return None

    @property
    def coeff_genset(self):
        return utils.poly_ring(utils.NoneVar)

    @cache
    def cached_product(self, u, v, basis2):
        if self == basis2:
            initial_dict = py.schubmult_q_fast({u: S.One}, v)
        else:
            initial_dict = {k: xreplace_genvars(x, self.coeff_genset, basis2.coeff_genset) for k, x in yz.schubmult_q_double_pair_generic(u, v).items()}
        return self.process_coeff_dict(initial_dict)

    @cache
    def cached_positive_product(self, u, v, basis2):
        return self.cached_product(u, v, basis2)

    def __call__(self, x):
        genset = self.genset
        # logger.debug(f"{x=} {type(x)=}")
        if not genset:
            genset = self.genset
        if isinstance(x, list) or isinstance(x, tuple):
            perm = Permutation(x)
            if not is_parabolic(perm, self.parabolic_index):
                raise ValueError(f"Permutation must be parabolic: {perm} is not")
            elem = self._from_dict({perm: 1})
        elif isinstance(x, Permutation):
            if not is_parabolic(x, self.parabolic_index):
                raise ValueError(f"Permutation must be parabolic: {x} is not")
            elem = self._from_dict({x: 1})
        elif isinstance(x, ParabolicQuantumDoubleSchubertAlgebraElement):
            return x
        else:
            dct = self.classical_basis(x)
            elem = 0
            if not isinstance(dct, spr.BasisSchubertAlgebraElement):
                return dct
            try:
                for k, v in dct.coeff_dict.items():
                    if elem == 0:
                        elem = v * self.classical_in_basis(k)
                    else:
                        elem += v * self.classical_in_basis(k)
            except ValueError:
                raise ValueError(f"Could not convert {x=} to quantum parabolic")
        return elem


def make_single_parabolic_quantum_basis(index_comp):
    return ParabolicQuantumSchubertAlgebraElement_basis(GeneratingSet("x"), index_comp)


@cache
def QPSx(*args):
    return make_single_parabolic_quantum_basis(args)
