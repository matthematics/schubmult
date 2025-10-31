from bisect import bisect_left
from functools import cache

import schubmult.mult.quantum as py
import schubmult.mult.quantum_double as yz
import schubmult.rings.schubert_ring as spr
import schubmult.utils.schub_lib as schub_lib
from schubmult.schub_lib.perm_lib import Permutation, longest_element, uncode
from schubmult.schub_lib.schub_poly import schubpoly_from_elems
from schubmult.symbolic import Add, Mul, Pow, S, Symbol, expand, expand_func, sympify
from schubmult.symmetric_polynomials import FactorialElemSym, QFactorialElemSym, coeffvars, degree, genvars, is_of_func_type, numvars
from schubmult.utils.logging import get_logger
from schubmult.utils.perm_utils import is_parabolic

from .abstract_schub_poly import PQDSchubPoly, QDSchubPoly
from .base_schubert_ring import BaseSchubertElement, BaseSchubertRing
from .poly_lib import complete_sym_poly, elem_sym_poly, elem_sym_poly_q, xreplace_genvars
from .variables import GeneratingSet, GeneratingSet_base, genset_dict_from_expr, poly_genset

q_var = GeneratingSet("q")

logger = get_logger(__name__)


def is_fact_elem_sym(obj):
    return is_of_func_type(obj, QFactorialElemSym)


class QuantumDoubleSchubertElement(BaseSchubertElement):
    def subs(self, old, new):
        return self.as_classical().subs(old, new).as_quantum()


class ParabolicQuantumDoubleSchubertElement(BaseSchubertElement):
    @property
    def index_comp(self):
        return self.ring.index_comp

    def kill_ideal(self):
        length = sum(self.index_comp)
        new_dict = {}
        for k, v in self.items():
            if len(k) <= length:
                new_dict[k] = v
        return self.ring.from_dict(new_dict)


class QuantumDoubleSchubertRing(BaseSchubertRing):
    def __hash__(self):
        return hash((self.genset, self.coeff_genset, "QDBS"))

    def __init__(self, genset, coeff_genset):
        super().__init__(genset, coeff_genset)
        self.dtype = type("QuantumDoubleSchubertElement", (QuantumDoubleSchubertElement,), {"ring": self})

    def __str__(self):
        return f"Quantum Double Schubert polynomial ring in {self.genset.label} and {self.coeff_genset.label}"

    def printing_term(self, k):
        return QDSchubPoly(k, self.genset.label, self.coeff_genset.label)

    def _coerce_mul(self, other):
        if isinstance(other, BaseSchubertElement):
            if other.ring == self:
                return other
            if isinstance(other.ring, QuantumDoubleSchubertRing) and other.ring.genset == self.genset:
                return other
            if type(other.ring) is QuantumSingleSchubertRing:
                if self.genset == other.ring.genset:
                    newbasis = QuantumDoubleSchubertRing(self.genset, poly_genset(0))
                    return newbasis.from_dict(other)
        return None

    def _coerce_add(self, other):
        if isinstance(other, BaseSchubertElement):
            if type(other.ring) is type(self):
                if self.genset == other.ring.genset and self.coeff_genset == other.ring.coeff_genset:
                    return other
        return None

    @property
    def symbol_elem_func(self):
        return QFactorialElemSym

    def elem_sym_subs(self, kk):
        elems = []
        for k in range(1, kk + 1):
            for p in range(1, k + 1):
                elems += [(Symbol(f"e_{p}_{k}"), elem_sym_poly_q(p, k, self.genset[1:], poly_genset(0)))]
        return dict(elems)

    @property
    def elem_sym(self):
        return QFactorialElemSym

    def is_elem_mul_type(self, other):
        return is_fact_elem_sym(other)

    # specifically symengine
    def elem_mul(self, ring_elem, elem):
        elem = sympify(elem)
        indexes = [self.genset.index(a) for a in genvars(elem)]
        ret = self.zero
        for k, v in ring_elem.items():
            perm_list = schub_lib.elem_sym_positional_perms_q(k, degree(elem), *indexes, q_var=q_var)
            # print(perm_list)
            for perm, df, sign, mul_val in perm_list:
                remaining_vars = [self.coeff_genset[perm[i - 1]] for i in indexes if perm[i - 1] == k[i - 1]]
                coeff = FactorialElemSym(degree(elem) - df, numvars(elem) - df, remaining_vars, coeffvars(elem))  # we need FactorialElemSym here
                ret += (v * sign * expand_func(coeff) * mul_val) * self(perm)
        return ret

    @cache
    def cached_product(self, u, v, basis2):
        return {k: xreplace_genvars(x, self.coeff_genset, basis2.coeff_genset) for k, x in yz.schubmult_q_double_pair_generic(u, v).items()}

    def in_quantum_basis(self, elem):
        return elem

    def in_classical_basis(self, elem):
        result = None
        for k, v in elem.items():
            if not result:
                result = v * self.quantum_as_classical_schubpoly(k)
            else:
                result += v * self.quantum_as_classical_schubpoly(k)
        return result if result else self.zero

    @property
    def classical_elem_func(self):
        basis = spr.DoubleSchubertRing(self.genset, self.coeff_genset)
        q_var = yz._vars.q_var

        def elem_func(p, k, varl1, varl2):
            if p == 0 and k >= 0:
                return basis.one
            if p < 0 or p > k:
                return basis.zero
            return (varl1[k - 1] - varl2[k - p]) * elem_func(p - 1, k - 1, varl1, varl2) + elem_func(p, k - 1, varl1, varl2) + q_var[k - 1] * elem_func(p - 2, k - 2, varl1, varl2)

        return elem_func

    @cache
    def quantum_as_classical_schubpoly(self, perm):
        return schubpoly_from_elems(perm, self.genset, self.coeff_genset, self.classical_elem_func)

    @cache
    def cached_schubpoly(self, k):
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

    def positive_elem_sym_rep(self, perm, index=1):
        if perm.inv == 0:
            return S.One
        ret = S.Zero
        L = schub_lib.pull_out_var(1, ~perm)
        for index_list, new_perm in L:
            ret += self.elem_sym(len(index_list), len(index_list), [self.genset[index2] for index2 in index_list], [self.coeff_genset[index]]) * self.positive_elem_sym_rep(~new_perm, index + 1)
        return ret

    def positive_elem_sym_rep_backward(self, perm):
        if perm.inv == 0:
            return S.One
        ret = S.Zero
        index = max((~perm).descents()) + 1
        L = schub_lib.pull_out_var(index, ~perm)
        for index_list, new_perm in L:
            ret += self.positive_elem_sym_rep_backward(~new_perm) * self.elem_sym(len(index_list), len(index_list), [self.genset[index2] for index2 in index_list], [self.coeff_genset[index]])
        return ret

    @property
    def mult_poly_double(self):
        return yz.mult_poly_q_double

    def from_expr(self, expr):
        ret = self.zero
        try:
            expr = sympify(expr)
            while expr != S.Zero:
                dct = genset_dict_from_expr(expr, self.genset)
                key = sorted(dct.keys(), reverse=True)[0]
                term = self.from_dict({uncode(key): dct[key]})
                ret += term
                expr -= term.as_polynomial()
                expr = expand(expr)
            return ret
        except Exception:
            return self.mul_expr(self.one, expr)

    def handle_sympoly(self, other):
        return expand_func(other)

    def mul_expr(self, elem, x):
        x = sympify(x)
        _Add = Add
        _Mul = Mul
        _Pow = Pow
        ind = self.genset.index(x)
        if ind != -1:
            return self.from_dict(yz.mult_poly_q_double(elem, x, self.genset, self.coeff_genset))
        if is_fact_elem_sym(x):
            # print(f"moo {x=}")
            if all(self.genset.index(a) != -1 for a in genvars(x)) and not any(self.genset.index(a) != -1 for a in coeffvars(x)):
                return self.elem_mul(elem, x)

            gens_to_remove = [a for a in genvars(x) if a not in self.genset]
            if any(self.genset.index(a) != -1 for a in genvars(x)) and len(gens_to_remove):
                return self.mul_expr(elem, x.split_out_vars(gens_to_remove))

            coeffs_to_remove = [a for a in coeffvars(x) if a in self.genset]

            if any(a in self.genset for a in coeffvars(x)) and len(coeffs_to_remove):
                return self.mul_expr(elem.split_out_vars(x.to_complete_sym(), coeffs_to_remove))
            return self.from_dict({k: (self.handle_sympoly(x)) * v for k, v in elem.items()})
        if isinstance(x, _Add):
            return self.sum([self.mul_expr(elem, arg) for arg in x.args])
        if isinstance(x, _Mul):
            res = elem
            for arg in x.args:
                res = self.mul_expr(res, arg)
            return res
        if isinstance(x, _Pow):
            res = elem
            for _ in range(int(x.args[1])):
                res = self.mul_expr(res, x.args[0])
            return res
        return self.from_dict({k: v * self.domain_new(x) for k, v in elem.items()})

    def new(self, x):
        genset = self.genset
        if not genset:
            genset = self.genset
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        if isinstance(x, list) or isinstance(x, tuple):
            elem = self.from_dict({Permutation(x): self.domain.one})
        elif isinstance(x, Permutation):
            elem = self.from_dict({x: self.domain.one})
        else:
            elem = self.from_expr(x)
        return elem


def QDSx(x, genset=GeneratingSet("y")):
    if isinstance(genset, str):
        genset = GeneratingSet(genset)
    return QuantumDoubleSchubertRing(GeneratingSet("x"), genset)(x)


class QuantumSingleSchubertRing(QuantumDoubleSchubertRing):
    def __init__(cls, genset):
        super().__init__(genset, poly_genset(0))

    def __hash__(self):
        return hash((self.genset, self.coeff_genset, "QBS"))

    def quantize(self, poly):
        r = spr.SingleSchubertRing(self.genset)
        ssp = r.from_expr(poly)
        sspq = self.from_dict(dict(ssp.items()))
        return sspq.as_polynomial()

    def _coerce_mul(self, other):
        """Coerce a basis schubert algebra element so it can be multiplied

        Args:
            other (_type_): _description_

        Returns:
            _type_: _description_
        """
        if other.ring == self:
            return other
        if type(other.ring) is QuantumDoubleSchubertRing:
            if self.genset == other.ring.genset:
                return other
        return None

    @cache
    def cached_product(self, u, v, basis2):
        if self == basis2:
            return py.schubmult_q_fast({u: S.One}, v)
        return {k: xreplace_genvars(x, self.coeff_genset, basis2.coeff_genset) for k, x in yz.schubmult_q_double_pair_generic(u, v).items()}

    @cache
    def cached_positive_product(self, u, v, basis2):
        return self.cached_product(u, v, basis2)

    def mul_expr(self, elem, x):
        x = sympify(x)
        _Add = Add
        _Mul = Mul
        _Pow = Pow
        ind = self.genset.index(x)
        if ind != -1:
            return self.from_dict(py.mult_poly_q(elem, x, self.genset))
        if isinstance(x, _Add):
            return self.sum([self.mul_expr(elem, arg) for arg in x.args])
        if isinstance(x, _Mul):
            res = elem
            for arg in x.args:
                res = self.mul_expr(res, arg)
            return res
        if isinstance(x, _Pow):
            res = elem
            for _ in range(int(x.args[1])):
                res = self.mul_expr(res, x.args[0])
            return res
        return self.from_dict({k: v * self.domain_new(x) for k, v in elem.items()})

    # def from_expr(self, x):
    #     x = sympify(x)
    #     result = py.mult_poly_q({Permutation([]): 1}, x, self.genset)
    #     return self.from_dict(result)

    def new(self, x):
        genset = self.genset
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        if isinstance(x, list) or isinstance(x, tuple):
            elem = self.from_dict({Permutation(x): self.domain.one})
        elif isinstance(x, Permutation):
            elem = self.from_dict({x: self.domain.one})
        elif isinstance(x, spr.DoubleSchubertElement):
            if x.genset == self.genset:
                return x.as_quantum()
        elif isinstance(x, ParabolicQuantumDoubleSchubertElement):
            return x.as_quantum()
        else:
            elem = self.from_expr(x)
        return elem

    # def in_SEM_basis(self):
    #     result = S.Zero
    #     for k, v in self.items():
    #         if len(k) > len(self._longest):
    #             parabolic_index = []
    #             start = 0
    #             index_comp = [*self._n, len(k) + 1 - self._N[-1]]
    #             for i in range(len(index_comp)):
    #                 end = start + index_comp[i]
    #                 parabolic_index += list(range(start + 1, end))
    #                 start = end
    #             otherlong = Permutation(list(range(parabolic_index[-1] + 1, 0, -1)))
    #             longpar = Permutation(longest_element(parabolic_index))
    #             longest = otherlong * longpar
    #         else:
    #             longest = self._longest
    #         result += v * schubpoly_from_elems(k, self.genset, self.coeff_genset, elem_func=self.ring.symbol_elem_func, mumu=~longest)
    #     return result


class ParabolicQuantumDoubleSchubertRing(BaseSchubertRing):
    def __hash__(self):
        return hash((self.genset, self.coeff_genset, self.index_comp, "PBGBG"))

    def __init__(self, genset, coeff_genset, index_comp):
        super().__init__(genset, coeff_genset)
        self._quantum_basis = QuantumDoubleSchubertRing(genset, coeff_genset)
        self._classical_basis = spr.DoubleSchubertRing(genset, coeff_genset)
        self._index_comp = index_comp
        self._n = list(index_comp)
        self._N = [sum(self._n[:i]) for i in range(len(self._n) + 1)]
        parabolic_index = []
        start = 0
        for i in range(len(index_comp)):
            end = start + index_comp[i]
            parabolic_index += list(range(start + 1, end))
            start = end
        self._parabolic_index = parabolic_index
        self._otherlong = Permutation(list(range(self._N[-1], 0, -1)))
        self._longest = self._otherlong * longest_element(parabolic_index)
        self.dtype = type("ParabolicQuantumDoubleSchubertElement", (ParabolicQuantumDoubleSchubertElement,), {"ring": self})

    def _coerce_mul(self, other):
        if isinstance(other, BaseSchubertElement):
            if other.ring == self:
                return other
            if isinstance(other.ring, ParabolicQuantumDoubleSchubertRing) and other.ring.genset == self.genset and self.index_comp == other.ring.index_comp:
                return other
            if type(other.ring) is ParabolicQuantumSingleSchubertRing:
                if self.genset == other.ring.genset and self.index_comp == other.ring.index_comp:
                    newbasis = ParabolicQuantumDoubleSchubertRing(self.genset, poly_genset(0), self.index_comp)
                    return newbasis.from_dict(other)
        return None

    def _coerce_add(self, other):
        if isinstance(other, BaseSchubertElement):
            if type(other.ring) is type(self):
                if self.genset == other.ring.genset and self.coeff_genset == other.ring.coeff_genset and self.index_comp == other.ring.index_comp:
                    return other
        return None

    @property
    def symbol_elem_func(self):
        def elem_func(p, k, varl1, varl2):  # noqa: ARG001
            if p == 0 and k >= 0:
                return 1
            if p < 0 or p > k:
                return 0
            return Add(*[(Symbol(f"e_{p - i}_{k}") if p - i > 0 else 1) * complete_sym_poly(i, k + 1 - p, [-v for v in varl2]) for i in range(p + 1)])

        return elem_func

    def elem_sym_subs(self, kk):
        elems = []
        elem_func = self.elem_sym
        for k in range(1, kk + 1):
            for p in range(1, k + 1):
                elems += [(Symbol(f"e_{p}_{k}"), elem_func(p, k, self.genset[1:], poly_genset(0)))]
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

    def elem_sym(self, p, k, varl1, varl2):
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

    @property
    def index_comp(self):
        return self._index_comp

    def process_coeff_dict(self, coeff_dict):
        max_len = max(len(w) for w in coeff_dict)
        parabolic_index = [*self._parabolic_index]
        if max_len > len(self._longest):
            parabolic_index = []
            start = 0
            index_comp = [*self._n, max_len + 1 - self._N[-1]]
            for i in range(len(index_comp)):
                end = start + index_comp[i]
                parabolic_index += list(range(start + 1, end))
                start = end
        return yz.apply_peterson_woodward(coeff_dict, parabolic_index)

    @cache
    def cached_product(self, u, v, basis2):
        initial_dict = {k: xreplace_genvars(x, self.coeff_genset, basis2.coeff_genset) for k, x in yz.schubmult_q_double_pair_generic(u, v).items()}
        return self.process_coeff_dict(initial_dict)

    def in_quantum_basis(self, elem):
        result = S.Zero
        for k, v in elem.items():
            result += v * schubpoly_from_elems(k, self.genset, self.coeff_genset, self.quantum_elem_func)
        return result

    def in_classical_basis(self, elem):
        result = S.Zero
        for k, v in elem.items():
            result += v * self.quantum_as_classical_schubpoly(k)
        return result

    # TODO: speed this up
    @cache
    def classical_in_basis(self, k):
        a = self.classical_basis(k)
        b = self(k)
        if expand(a.as_polynomial() - b.as_polynomial()) == S.Zero:
            return b

        cd = b.as_classical()
        for k2, v in cd.items():
            if expand(v) == S.Zero:
                continue
            if k != k2:
                b -= v * self.classical_in_basis(k2)
        return self.from_dict({k: self.domain_new(v) for k, v in b.items() if expand(v) != S.Zero})

    @property
    def classical_elem_func(self):
        basis = spr.DoubleSchubertRing(self.genset, self.coeff_genset)
        q_var = yz._vars.q_var

        def elem_func(p, k, varl1, varl2):
            if p == 0 and k >= 0:
                return basis.one
            if p < 0 or p > k:
                return basis.zero
            if k <= self._N[1]:
                return basis.elem_sym(p, k, varl1, varl2)  # check: this may speed it up
            ret = basis.zero
            j = bisect_left(self._N, k)
            if j < len(self._N) and k == self._N[j]:
                ret = (-((S.NegativeOne) ** (self._n[j - 1]))) * q_var[j - 1] * elem_func(p - self._N[j] + self._N[j - 2], self._N[j - 2], varl1, varl2)
            ret += elem_func(p, k - 1, varl1, varl2) + basis.elem_sym(1, 1, varl1[k - 1], varl2[k - p]) * elem_func(p - 1, k - 1, varl1, varl2)
            return ret

        return elem_func

    @property
    def quantum_elem_func(self):
        basis = QuantumDoubleSchubertRing(self.genset, self.coeff_genset)
        q_var = yz._vars.q_var

        def elem_func(p, k, varl1, varl2):
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
            return ret

        return elem_func

    def printing_term(self, k):
        return PQDSchubPoly(k, self.genset.label, self.coeff_genset.label, self.index_comp)

    @cache
    def quantum_as_classical_schubpoly(self, perm):
        k = perm
        if len(k) > len(self._longest):
            parabolic_index = []
            start = 0
            index_comp = [*self._n, len(k) + 1 - self._N[-1]]
            for i in range(len(index_comp)):
                end = start + index_comp[i]
                parabolic_index += list(range(start + 1, end))
                start = end
            otherlong = Permutation(list(range(parabolic_index[-1] + 1, 0, -1)))
            longpar = Permutation(longest_element(parabolic_index))
            longest = otherlong * longpar
        else:
            longest = self._longest
        return schubpoly_from_elems(perm, self.genset, self.coeff_genset, elem_func=self.classical_elem_func, mumu=~longest)

    @cache
    def cached_schubpoly(self, k):
        if len(k) > len(self._longest):
            parabolic_index = []
            start = 0
            index_comp = [*self._n, len(k) + 1 - self._N[-1]]
            for i in range(len(index_comp)):
                end = start + index_comp[i]
                parabolic_index += list(range(start + 1, end))
                start = end
            otherlong = Permutation(list(range(parabolic_index[-1] + 1, 0, -1)))
            longpar = Permutation(longest_element(parabolic_index))
            longest = otherlong * longpar
        else:
            longest = self._longest
        return schubpoly_from_elems(k, self.genset, self.coeff_genset, elem_func=self.elem_sym, mumu=~longest)

    @cache
    def cached_positive_product(self, u, v, basis2):
        initial_dict = {k: xreplace_genvars(x, self.coeff_genset, basis2.coeff_genset if basis2.coeff_genset else poly_genset(0)) for k, x in yz.schubmult_q_generic_partial_posify(u, v).items()}
        return self.process_coeff_dict(initial_dict)

    @property
    def double_mul(self):
        from schubmult.mult.quantum_double import _vars

        def do_double_mul(perm_dict, v, var2=_vars.var2, var3=_vars.var3, q_var=_vars.q_var):
            coeff_dict = yz.schubmult_q_double_fast(perm_dict, v, var2, var3, q_var)
            return self.process_coeff_dict(coeff_dict)

        return do_double_mul

    @property
    def single_mul(self):
        from schubmult.mult.quantum_double import _vars

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

    def from_expr(self, expr):
        ret = self.zero
        try:
            expr = sympify(expr)
            while expr != S.Zero:
                dct = genset_dict_from_expr(expr, self.genset)
                key = sorted(dct.keys(), reverse=True)[0]
                new_key = []
                coeff1 = dct[key]
                for i in range(len(self.index_comp)):
                    new_key += sorted(key[len(new_key) : len(new_key) + self.index_comp[i]])
                new_key = tuple(new_key)
                try:
                    coeff2 = dct[new_key]
                except KeyError:
                    raise ValueError(f"{expr} does not have symmetry properties")
                if expand(coeff1 - coeff2) != S.Zero:
                    raise ValueError(f"{expr} does not have symmetry properties")
                term = self.from_dict({uncode(new_key): dct[new_key]})
                ret += term
                expr -= term.as_polynomial()
                expr = expand(expr)
            return ret
        except Exception:
            raise

    # TODO: speed this up
    def mul_expr(self, elem, x):
        # dct = self.classical_basis.mul_expr(self.in_classical_basis(elem), x)
        # elem = self.zero
        # if not isinstance(dct, BaseSchubertElement):
        #     return dct
        # for k, v in dct.items():
        #     if expand(v) == S.Zero:
        #         continue
        #     if elem == self.zero:
        #         elem = v * self.classical_in_basis(k)
        #     else:
        #         elem += v * self.classical_in_basis(k)
        return elem * self.from_expr(x)

    def __call__(self, x):
        genset = self.genset
        if not genset:
            genset = self.genset
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        if isinstance(x, list) or isinstance(x, tuple):
            perm = Permutation(x)
            if not is_parabolic(perm, self.parabolic_index):
                raise ValueError(f"Permutation must be parabolic: {perm} is not")
            elem = self.from_dict({perm: self.domain.one})
        elif isinstance(x, Permutation):
            if not is_parabolic(x, self.parabolic_index):
                raise ValueError(f"Permutation must be parabolic: {x} is not")
            elem = self.from_dict({x: self.domain.one})
        elif isinstance(x, ParabolicQuantumDoubleSchubertElement):
            return x
        else:
            elem = self.from_expr(x)
        return elem


QSx = QuantumSingleSchubertRing(GeneratingSet("x"))

QuantumDoubleSchubertPolynomial = QuantumDoubleSchubertElement


def make_parabolic_quantum_basis(index_comp, coeff_genset):
    return ParabolicQuantumDoubleSchubertRing(GeneratingSet("x"), coeff_genset, index_comp)


def QPDSx_index(*args):
    def this_QPDSx(x, coeff_genset=GeneratingSet("y")):
        return make_parabolic_quantum_basis(args, poly_genset(coeff_genset) if isinstance(coeff_genset, str) else coeff_genset)(x)

    return this_QPDSx


@cache
def QPDSx(*args):
    return QPDSx_index(*args)


class ParabolicQuantumSingleSchubertRing(ParabolicQuantumDoubleSchubertRing):
    def __init__(self, genset, index_comp):
        super().__init__(genset, poly_genset(0), index_comp)

    def __hash__(self):
        return hash((self.genset, self.coeff_genset, self.index_comp, "PQSB"))

    @cache
    def cached_schubpoly(self, k):
        if len(k) > len(self._longest):
            parabolic_index = []
            start = 0
            index_comp = [*self._n, len(k) + 1 - self._N[-1]]
            for i in range(len(index_comp)):
                end = start + index_comp[i]
                parabolic_index += list(range(start + 1, end))
                start = end
            otherlong = Permutation(list(range(parabolic_index[-1] + 1, 0, -1)))
            longpar = Permutation(longest_element(parabolic_index))
            longest = otherlong * longpar
        else:
            longest = self._longest
        return schubpoly_from_elems(k, self.genset, self.coeff_genset, elem_func=self.elem_sym, mumu=~longest)

    def elem_sym(self, p, k, varl1, varl2):
        if p < 0 or p > k:
            return 0
        if p == 0 and k >= 0:
            return 1
        if k <= self._N[1]:
            return elem_sym_poly(p, k, varl1, poly_genset(0))
        ret = 0
        j = bisect_left(self._N, k)
        if j < len(self._N) and k == self._N[j]:
            ret = (-((-1) ** (self._n[j - 1]))) * q_var[j - 1] * self.elem_sym(p - self._N[j] + self._N[j - 2], self._N[j - 2], varl1, varl2)
        ret += self.elem_sym(p, k - 1, varl1, varl2) + varl1[k - 1] * self.elem_sym(p - 1, k - 1, varl1, varl2)
        return ret

    def _coerce_mul(self, other):
        if isinstance(other, BaseSchubertElement):
            if other.ring == self:
                if self.genset == other.ring.genset:
                    return other
            if type(other.ring) is ParabolicQuantumDoubleSchubertRing:
                if self.genset == other.ring.genset and self.index_comp == other.ring.index_comp:
                    newbasis = ParabolicQuantumDoubleSchubertRing(self.genset, poly_genset(0), self.index_comp)
                    return newbasis.from_dict(other)
        return None

    @property
    def coeff_genset(self):
        return poly_genset(0)

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
        if not genset:
            genset = self.genset
        if isinstance(x, list) or isinstance(x, tuple):
            perm = Permutation(x)
            if not is_parabolic(perm, self.parabolic_index):
                raise ValueError(f"Permutation must be parabolic: {perm} is not")
            elem = self.from_dict({perm: self.domain.one})
        elif isinstance(x, Permutation):
            if not is_parabolic(x, self.parabolic_index):
                raise ValueError(f"Permutation must be parabolic: {x} is not")
            elem = self.from_dict({x: self.domain.one})
        elif isinstance(x, ParabolicQuantumDoubleSchubertElement):
            return x
        else:
            elem = self.from_expr(x)
        return elem


def make_single_parabolic_quantum_basis(index_comp):
    return ParabolicQuantumSingleSchubertRing(GeneratingSet("x"), index_comp)


@cache
def QPSx(*args):
    return make_single_parabolic_quantum_basis(args)
