from bisect import bisect_left
from functools import cache

import schubmult.mult.quantum as py
import schubmult.mult.quantum_double as yz
import schubmult.rings.schubert.schubert_ring as spr
from schubmult.combinatorial_reps.permutation import Permutation, uncode
from schubmult.symbolic import Add, S, Symbol, expand, sympify
from schubmult.symbolic.poly.poly_lib import complete_sym_poly, elem_sym_poly, xreplace_genvars
from schubmult.symbolic.poly.schub_poly import schubpoly_from_elems
from schubmult.symbolic.poly.variables import GeneratingSet, GeneratingSet_base, genset_dict_from_expr, poly_genset
from schubmult.utils.perm_utils import is_parabolic

from ..printing import PQDSchubPoly
from .base_schubert_ring import BaseSchubertElement, BaseSchubertRing
from .quantum_double_schubert_ring import QuantumDoubleSchubertRing

q_var = GeneratingSet("q")


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
        self._longest = self._otherlong * Permutation.longest_element(*parabolic_index)
        self.dtype = type("ParabolicQuantumDoubleSchubertElement", (ParabolicQuantumDoubleSchubertElement,), {"ring": self})

    def _coerce_mul(self, other):
        from .parabolic_quantum_schubert_ring import ParabolicQuantumSingleSchubertRing

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
        qv = yz._vars.q_var

        def elem_func(p, k, varl1, varl2):
            if p == 0 and k >= 0:
                return basis.one
            if p < 0 or p > k:
                return basis.zero
            if k <= self._N[1]:
                return basis.elem_sym(p, k, varl1, varl2)
            ret = basis.zero
            j = bisect_left(self._N, k)
            if j < len(self._N) and k == self._N[j]:
                ret = (-((S.NegativeOne) ** (self._n[j - 1]))) * qv[j - 1] * elem_func(p - self._N[j] + self._N[j - 2], self._N[j - 2], varl1, varl2)
            ret += elem_func(p, k - 1, varl1, varl2) + basis.elem_sym(1, 1, varl1[k - 1], varl2[k - p]) * elem_func(p - 1, k - 1, varl1, varl2)
            return ret

        return elem_func

    @property
    def quantum_elem_func(self):
        basis = QuantumDoubleSchubertRing(self.genset, self.coeff_genset)
        qv = yz._vars.q_var

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
                ret = (-((-1) ** (self._n[j - 1]))) * qv[j - 1] * elem_func(p - self._N[j] + self._N[j - 2], self._N[j - 2], varl1, varl2)
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
            longpar = Permutation.longest_element(*parabolic_index)
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
            longpar = Permutation.longest_element(*parabolic_index)
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

    def mul_expr(self, elem, x):
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


def make_parabolic_quantum_basis(index_comp, coeff_genset):
    return ParabolicQuantumDoubleSchubertRing(GeneratingSet("x"), coeff_genset, index_comp)


def QPDSx_index(*args):
    def this_QPDSx(x, coeff_genset=GeneratingSet("y")):
        return make_parabolic_quantum_basis(args, poly_genset(coeff_genset) if isinstance(coeff_genset, str) else coeff_genset)(x)

    return this_QPDSx


@cache
def QPDSx(*args):
    return QPDSx_index(*args)
