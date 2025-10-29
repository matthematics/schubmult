from functools import cache, cached_property

import schubmult.mult.double as yz
import schubmult.mult.positivity as pos
import schubmult.mult.single as py
import schubmult.rings.abstract_schub_poly as spolymod
import schubmult.rings.quantum_schubert_ring as qsr
import schubmult.utils.schub_lib as schub_lib
from schubmult.schub_lib.perm_lib import Permutation, uncode
from schubmult.schub_lib.schub_poly import schubpoly_classical_from_elems, schubpoly_from_elems
from schubmult.symbolic import Add, DomainElement, Mul, Pow, S, Symbol, expand, expand_func, is_of_func_type, sympify, sympify_sympy
from schubmult.symmetric_polynomials import CompleteSym_base, ElemSym, ElemSym_base, FactorialElemSym, coeffvars, degree, genvars, numvars, split_out_vars
from schubmult.utils.logging import get_logger
from schubmult.utils.perm_utils import add_perm_dict

from .base_schubert_ring import BaseSchubertElement, BaseSchubertRing
from .poly_lib import elem_sym_poly, xreplace_genvars
from .tensor_ring import TensorRing
from .variables import CustomGeneratingSet, GeneratingSet, GeneratingSet_base, MaskedGeneratingSet, NotEnoughGeneratorsError, genset_dict_from_expr, poly_genset

logger = get_logger(__name__)


def is_fact_elem_sym(obj):
    return is_of_func_type(obj, ElemSym_base)


def is_fact_complete_sym(obj):
    return is_of_func_type(obj, CompleteSym_base)


class DoubleSchubertElement(BaseSchubertElement):
    """Algebra with sympy coefficients
    and a dict basis
    """

    def to_genset_dict(self, trim=False):
        gdict = genset_dict_from_expr(self.as_polynomial(), self.ring.genset)

        def _trim_tuple(tup):
            tup = [*tup]
            while len(tup) > 0 and tup[-1] == 0:
                tup.pop()
            return (*tup,)

        if trim:
            new_dict = {}
            for flop, val in gdict.items():
                real_tup = _trim_tuple(flop)
                new_dict[real_tup] = new_dict.get(real_tup, S.Zero) + val
        else:
            new_dict = gdict
        return new_dict

    def divdiff(self, i):
        return self.ring.from_dict({k.swap(i - 1, i): v for k, v in self.items() if i - 1 in k.descents()})

    def simpleref(self, i):
        return self + self.divdiff(i).mult_poly(self.ring.genset[i + 1] - self.ring.genset[i])

    def act(self, perm):
        perm = Permutation(perm)
        dset = perm.descents()
        if len(dset) == 0:
            return self
        i = next(iter(dset))
        return self.simpleref(i + 1).act(perm.swap(i, i + 1))

    def max_index(self):
        return max([max([0, *list(k.descents(zero_indexed=False))]) for k in self.keys()])

    def eval(self, x):
        ret = self
        for v, val in x.items():
            ret = ret.pull_out_gen(v)
            ret = ret.ring.from_dict({k: v2.subs(v, val) for k, v2 in ret.items()})
        if len(ret.keys()) == 1 and next(iter(ret.keys())) == Permutation([]):
            return ret[Permutation([])]
        return ret

    # def to_pos_elem_sym(self):
    #     ret = self
    #     for v, val in x.items():
    #         ret = ret.pull_out_gen(v)
    #         ret = ret.ring.from_dict({k: v2.subs(v, val) for k, v2 in ret.items()})
    #     if len(ret.keys()) == 1 and next(iter(ret.keys())) == Permutation([]):
    #         return ret[Permutation([])]
    #     return ret

    def subs(self, old, new):
        result = 0
        if self.ring.genset.index(old) != -1:
            result = 0
            index = self.ring.genset.index(old)
            mindex = self.max_index()
            if mindex < index:
                return self
            perm = Permutation([]).swap(index - 1, mindex)  # index to max index + 1
            transf = self.act(perm)
            for k, v in transf.items():
                perm = k
                coeff_gens = self.ring.coeff_genset
                L = schub_lib.pull_out_var(mindex + 1, perm)
                for index_list, new_perm in L:
                    result += self.ring.from_dict({new_perm: v}).mult_poly(Mul(*[(new - coeff_gens[index2]) for index2 in index_list]))
            return result

        for k, v in self.items():
            if self.ring.coeff_genset.label is None:
                add_dict = {k: v.subs(old, new)}
            else:
                coeff_genset = self.ring.coeff_genset
                if coeff_genset.index(old) != -1:
                    genset_list = [coeff_genset[i] for i in range(len(coeff_genset))]
                    genset_list[coeff_genset.index(old)] = 0
                    custom_genset = CustomGeneratingSet(genset_list)
                    new_add_dict = {k2: sympify(v2).subs(old, new) for k2, v2 in yz.schubmult_double({(): v}, k, custom_genset, coeff_genset).items()}  # remove the variable
                    add_dict = {}
                    for k3, v3 in new_add_dict.items():
                        to_add_dict = yz.schubmult_double({(): v3}, k3, coeff_genset, custom_genset)
                        add_dict = add_perm_dict(add_dict, to_add_dict)
                else:
                    add_dict = {k: sympify(v).subs(old, new)}
            for k5, v5 in add_dict.items():
                if any(self.ring.genset.index(s) != -1 for s in sympify(v5).free_symbols):
                    result += self.ring.from_dict({k5: 1}).mult_poly(v5)
                else:
                    result += self.ring.from_dict({k5: v5})
        return result

    @property
    def free_symbols(self):
        ret = set()
        for k, v in self.items():
            ret.update(v.free_symbols)
            perm = k
            if len(perm.descents()) > 0:
                ret.update([self.ring.genset[i] for i in range(1, max(perm.descents()) + 2)])
            if self.ring.coeff_genset:
                genset2 = self.ring.coeff_genset
                perm2 = ~perm
                if len(perm2.descents()) > 0:
                    ret.update([genset2[i] for i in range(1, max(perm2.descents()) + 2)])
        return ret

    def pull_out_gen(self, gen):
        ind = self.ring.genset.index(gen)
        if ind == -1:
            ind = self.ring.coeff_genset.index(gen)
            if ind == -1:
                raise ValueError(f"{gen} passed but is not a generator")
            gens2 = MaskedGeneratingSet(self.ring.coeff_genset, [ind])
            gens2.set_label(f"({self.ring.coeff_genset.label}\\{gen})")
            new_basis = DoubleSchubertRing(self.ring.genset, gens2)
            ret = new_basis.zero
            for perm, val in self.items():
                L = schub_lib.pull_out_var(ind, ~perm)
                for index_list, new_perm in L:
                    toadd = S.One
                    # for index2 in index_list:
                    #     toadd *= self.ring.genset[index2] - gen
                    ret += FactorialElemSym(len(index_list), len(index_list), [self.ring.genset[index2] for index2 in index_list], [gen]) * val * new_basis(~new_perm)
            return ret
        gens2 = MaskedGeneratingSet(self.ring.genset, [ind])
        gens2.set_label(f"({self.ring.genset.label}\\{gen})")
        new_basis = DoubleSchubertRing(gens2, self.ring.coeff_genset)
        ret = new_basis.zero
        for perm, val in self.items():
            L = schub_lib.pull_out_var(ind, perm)
            for index_list, new_perm in L:
                toadd = S.One
                for index2 in index_list:
                    toadd *= gen - self.ring.coeff_genset[index2]
                ret += toadd * val * new_basis(new_perm)
        return ret

    def in_CEM_basis(self):
        result = S.Zero
        for k, v in self.items():
            result += sympify(v) * schubpoly_classical_from_elems(k, self.ring.genset, self.ring.coeff_genset, elem_func=self.ring.symbol_elem_func)
        return result

    def cem_rep(self, elem_func, mumu=None):
        result = S.Zero
        if mumu:
            for k, v in self.items():
                result += sympify(v) * schubpoly_from_elems(k, self.ring.genset, self.ring.coeff_genset, elem_func=elem_func, mumu=mumu)
        else:
            for k, v in self.items():
                result += sympify(v) * schubpoly_classical_from_elems(k, self.ring.genset, self.ring.coeff_genset, elem_func=elem_func)
        return result

    def coproduct(self, *indices, alt_coeff_genset=None, on_coeff_gens=False, gname1=None, gname2=None):
        result_dict = {}
        genset = self.ring.genset
        if on_coeff_gens:
            genset = self.ring.coeff_genset
        if gname1 is None:
            gname1 = f"{genset.label}_A"  # "("+", ".join([f"{genset.label}_{i}" for i in indices])+")"
        if gname2 is None:
            gname2 = f"{genset.label}_B"  # f"{genset.label}\\{{"+", ".join([f"{genset.label}_{i}" for i in indices])+"}"
        gens2 = MaskedGeneratingSet(genset, indices)
        gens1 = gens2.complement()
        gens1.set_label(gname1)
        gens2.set_label(gname2)
        for k, v in self.items():
            key = k
            if isinstance(self.ring, SingleSchubertRing) and not alt_coeff_genset:
                coprod_dict = py.schub_coprod_py(key, indices)
            else:
                if on_coeff_gens:
                    coprod_dict = yz.schub_coprod_double(~key, indices, self.ring.genset, alt_coeff_genset if alt_coeff_genset else self.ring.genset)
                else:
                    coprod_dict = yz.schub_coprod_double(key, indices, self.ring.coeff_genset, alt_coeff_genset if alt_coeff_genset else self.ring.coeff_genset)
            if on_coeff_gens:
                result_dict = add_perm_dict(result_dict, {(~k1, ~k2): v * v2 for (k1, k2), v2 in coprod_dict.items()})
            else:
                result_dict = add_perm_dict(result_dict, {k: v * v2 for k, v2 in coprod_dict.items()})
        if on_coeff_gens:
            basis = TensorRing(
                DoubleSchubertRing(self.ring.genset, gens1),
                DoubleSchubertRing(alt_coeff_genset if alt_coeff_genset else self.ring.genset, gens2),
            )
        else:
            basis = TensorRing(
                DoubleSchubertRing(gens1, self.ring.coeff_genset),
                DoubleSchubertRing(gens2, alt_coeff_genset if alt_coeff_genset else self.ring.coeff_genset),
            )
        return basis.from_dict(result_dict)

    @cached_property
    def max_gens(self):
        return max([max(k.descents()) for k in self.keys()])

    def positive_elem_sym_rep(self):
        res = S.Zero
        for k, val in self.items():
            res += val * self.ring.positive_elem_sym_rep(k)
        return res

    def positive_elem_sym_rep_backward(self):
        res = S.Zero
        for k, val in self.items():
            res += val * self.ring.positive_elem_sym_rep_backward(k)
        return res


# schubmult_double_down
class DoubleSchubertRing(BaseSchubertRing):
    def __hash__(self):
        return hash((self.genset, self.coeff_genset, "DBS"))

    def __init__(self, genset, coeff_genset, domain=None):
        super().__init__(genset, coeff_genset, domain)
        self.dtype = type("DoubleSchubertElement", (DoubleSchubertElement,), {"ring": self})

    def __str__(self):
        return f"Double Schubert polynomial ring in {self.genset.label} and {self.coeff_genset.label}"

    def rmul(self, elem, other):
        import schubmult.rings.free_algebra as fa
        import schubmult.rings.free_algebra_basis as fb

        if isinstance(other, fa.FreeAlgebraElement):
            other = other.change_basis(fb.SchubertBasis)
            manip = elem
            ring = self
            if not isinstance(self, SingleSchubertRing):
                manip = SingleSchubertRing(self.genset)([]) * elem
                ring = manip.ring
            ret0 = manip.ring.zero
            for (k, n), v in other.items():
                for k1, v1 in manip.items():
                    n2 = 0
                    if k1.inv > 0:
                        n2 = max(k1.descents(False))
                    if n > n2:
                        continue
                    toshift = n2 - n
                    new_k = uncode([*([0] * toshift), *k.code])
                    if (k1 * ~new_k).inv != k1.inv - new_k.inv:
                        continue
                    tosplit = ring(k1 * (~new_k))
                    dct = tosplit.coproduct(*list(range(1, toshift + 1)))
                    for (perm1, perm2), v2 in dct.items():
                        if perm2.inv == 0:
                            ret0 += v * v1 * v2 * ring(perm1)
            return self([]) * ret0
        try:
            other = self.domain_new(other)
            return self.from_dict({k: v * other for k, v in elem.items()})
        except Exception:
            return self.mul_expr(elem, other)

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

    # def mul(self, elem, other, _sympify=False):
    #     try:
    #         other = self.domain_new(other)
    #         return self.from_dict({k: other * v for k, v in elem.items()})
    #     except CoercionFailed:
    #         pass
    #     if isinstance(other, BaseSchubertElement):
    #         other = self._coerce_mul(other)
    #         if not other:
    #             raise CoercionFailed(f"Could not coerce {other} of type {type(other)} to {type(elem)}")
    #         return self.from_dict(utils._mul_schub_dicts(elem, other, elem.ring, other.ring, _sympify=_sympify))

    #     return self.mul_expr(elem, other)

    def printing_term(self, k, prefix=""):
        return spolymod.DSchubPoly(k, self.genset.label, self.coeff_genset.label, prefix=prefix)

    def _coerce_mul(self, other):
        if isinstance(other, BaseSchubertElement):
            if isinstance(other.ring, qsr.QuantumDoubleSchubertRing):
                return other.as_classical()
            if isinstance(other.ring, ElemDoubleSchubertRing):
                return other
            if isinstance(other.ring, DoubleSchubertRing):
                return other
        return None

    def _coerce_add(self, other):  # noqa: ARG002
        return None

    @property
    def elem_sym(self):
        return FactorialElemSym

    def is_elem_mul_type(self, other):
        return is_fact_elem_sym(other)

    # specifically symengine
    def elem_mul(self, ring_elem, elem):
        elem = sympify(elem)
        indexes = [self.genset.index(a) for a in genvars(elem)]
        ret = self.zero
        elem_sympy = sympify_sympy(elem)
        for k, v in ring_elem.items():
            perm_list = schub_lib.elem_sym_positional_perms(k, degree(elem), *indexes)
            for perm, df, sign in perm_list:
                remaining_vars = [self.coeff_genset[perm[i - 1]] for i in indexes if perm[i - 1] == k[i - 1]]
                coeff = elem_sympy.func(degree(elem) - df, numvars(elem) - df, remaining_vars, coeffvars(elem))  # leave as elem sym
                ret += (v * sign * expand_func(coeff)) * self(perm)
        return ret

    @property
    def symbol_elem_func(self):
        return FactorialElemSym

    def schubert_schur_elem_func(self, numvars):
        ring = self @ self

        def elem_func(p, k, *args):  # noqa: ARG001
            if p < 0:
                return ring.zero
            if p > k:
                return ring.zero
            if k >= 0 and p == 0:
                return ring.one
            if k >= numvars:
                return ring((uncode([0] * (k - p) + [1] * p), Permutation([])))
            return ring((Permutation([]), uncode([0] * (k - p) + [1] * p)))

        return elem_func

    def in_schubert_schur_basis(self, perm, numvars):
        elem_func = self.schubert_schur_elem_func(numvars)
        if perm.inv == 0:
            return elem_func(0, 0)
        extra = len(perm) - numvars
        dom = uncode([numvars] * extra + list(range(numvars - 1, 0, -1)))
        return schubpoly_from_elems(perm, self.genset, self.coeff_genset, elem_func=elem_func, mumu=dom)

    # def in_shifted_SEM_basis(self, shift):
    #     result = S.Zero
    #     for k, v in self.items():
    #         result += sympify(v) * schubpoly_from_elems(k, self.ring.genset, self.ring.coeff_genset, elem_func=self.ring.shifted_symbol_elem_func(shift))
    #     return result

    def in_descending_schur_basis(self, perm, numvars):
        if numvars == 1:
            if perm == uncode([1]) or perm.inv == 0:
                return Sx(perm)
            return self.zero
        mid_res = self.in_schubert_schur_basis(perm, numvars)
        new_ring = TensorRing(mid_res.ring, self)
        result = new_ring.zero
        for k, v in mid_res.items():
            second_part = self.in_descending_schur_basis(k[1], numvars - 1)
            if numvars == 2:
                for k1, v2 in second_part.items():
                    result += new_ring.from_dict({(k[0], k1): v * v2})
            else:
                for (k1, k2), v2 in second_part.items():
                    result += new_ring.from_dict({((k[0], k1), k2): v * v2})
        return result

    def elem_sym_subs(self, kk):
        elems = []
        for k in range(1, kk + 1):
            for p in range(1, k + 1):
                elems += [(Symbol(f"e_{p}_{k}"), elem_sym_poly(p, k, self.genset[1:], poly_genset(0)))]
        return dict(elems)

    @staticmethod
    def flip(elem):
        R = DoubleSchubertRing(CustomGeneratingSet([0, *coeffvars(elem)]), CustomGeneratingSet([0, *genvars(elem)]))
        p = degree(elem)
        K = numvars(elem) + 1 - p
        poly = R(uncode([*list((K - 1) * [0]), p]))
        # print(poly)
        return poly.in_CEM_basis()

    def in_quantum_basis(self, elem):
        result = S.Zero
        for k, v in elem.items():
            result += v * self.quantum_schubpoly(k)
        return result

    def in_classical_basis(self, elem):
        return elem

    @cache
    def quantum_schubpoly(self, perm):
        return schubpoly_classical_from_elems(perm, self.genset, self.coeff_genset, self.quantum_elem_func)

    @cache
    def cached_product(self, u, v, basis2):
        return {k: xreplace_genvars(x, self.coeff_genset, basis2.coeff_genset) for k, x in yz.schubmult_double_pair_generic(u, v).items()}

    @cache
    def cached_positive_product(self, u, v, basis2):
        return {k: xreplace_genvars(x, self.coeff_genset, basis2.coeff_genset) for k, x in pos.schubmult_generic_partial_posify(u, v).items()}

    @property
    def double_mul(self):
        return yz.schubmult_double

    @property
    def single_mul(self):
        return py.schubmult_py

    @property
    def mult_poly_single(self):
        return py.mult_poly_py

    @property
    def mult_poly_double(self):
        return yz.mult_poly_double

    @property
    def quantum_elem_func(self):
        basis = qsr.QuantumDoubleSchubertRing(self.genset, self.coeff_genset)

        def elem_func(p, k, varl1, varl2, xstart=0, ystart=0):
            if p > k:
                return basis(0)
            if p == 0:
                return basis([])
            if p == 1:
                res = basis(varl1[xstart] - varl2[ystart])
                for i in range(1, k):
                    res += basis(varl1[xstart + i] - varl2[ystart + i])
                return res
            if p == k:
                res = basis((varl1[xstart] - varl2[ystart]) * (varl1[xstart + 1] - varl2[ystart]))
                for i in range(2, k):
                    res *= basis(varl1[i + xstart] - varl2[ystart])
                return res
            mid = k // 2
            xsm = xstart + mid
            ysm = ystart + mid
            kmm = k - mid
            res = elem_func(p, mid, varl1, varl2, xstart, ystart) + elem_func(
                p,
                kmm,
                varl1,
                varl2,
                xsm,
                ysm,
            )
            for p2 in range(max(1, p - kmm), min(p, mid + 1)):
                res += elem_func(p2, mid, varl1, varl2, xstart, ystart) * elem_func(
                    p - p2,
                    kmm,
                    varl1,
                    varl2,
                    xsm,
                    ysm - p2,
                )
            return res

        return elem_func

    def monomial_schub(self, monom):
        monom = [*monom]
        while len(monom) > 0 and monom[-1] == 0:
            monom.pop()
        return self._monomial_schub_cache(tuple(monom))

    @cache
    def _monomial_schub_cache(self, monom):
        srt_perm = Permutation.sorting_perm([-i for i in monom])
        schub_perm = uncode(sorted(monom, reverse=True))
        return self.from_dict({(schub_perm, 0): S.One}).act(srt_perm)

    @cache
    def cached_schubpoly(self, k):
        return schubpoly_classical_from_elems(k, self.genset, self.coeff_genset, elem_func=elem_sym_poly)

    def complete_mul(self, elem, x):
        x = sympify(x)
        x_sympy = sympify_sympy(x)
        indexes = {self.genset.index(a) for a in genvars(x)}
        ret = self.zero
        for k, v in elem.items():
            perm_list = schub_lib.complete_sym_positional_perms(k, degree(x), *indexes)
            for perm, df, sign in perm_list:
                remaining_vars = [self.coeff_genset[perm[i - 1]] for i in {*indexes, *[j + 1 for j in range(len(perm)) if perm[j] != k[j]]}]
                coeff = x_sympy.func(degree(x) - df, numvars(x) + df, remaining_vars, coeffvars(x))  # leave as elem sym
                ret += (sign * v * expand_func(coeff)) * self(perm)
        return ret

    def handle_sympoly(self, other):
        return expand_func(other)

    def single_variable(self, elem, varnum):
        ret = self.zero
        for u, v in elem.items():
            ret += v * self.coeff_genset[u[varnum - 1]] * self(u)
            new_perms = schub_lib.elem_sym_positional_perms(u, 1, varnum)
            for perm, udiff, sign in new_perms:
                if udiff == 1:
                    ret += (sign * v) * self(perm)
        return ret

    def from_expr(self, expr):
        # ret = self.zero
        # try:
        #     expr = sympify(expr)
        #     while expr != S.Zero:
        #         dct = genset_dict_from_expr(expr, self.genset)
        #         key = sorted(dct.keys(), reverse=True)[0]
        #         term = self.from_dict({uncode(key): dct[key]})
        #         ret += term
        #         expr -= term.as_polynomial()
        #         expr = expand(expr)
        #     return ret
        # except Exception:
        return super().from_expr(expr)

    # def from_expr(self, expr):
    #     from .polynomial_algebra import PA
    #     ret = self.zero
    #     dexpr = PA.from_expr(expr)
    #     try:
    #         while dexpr.keys():
    #             # dct = genset_dict_from_expr(expr, self.genset)
    #             key = next(iter(sorted(dexpr.keys(), reverse=True)))
    #             # print(f"{key=} {dexpr[key]=}")
    #             term = self.from_dict({uncode(key): dexpr[key]})
    #             # print(f"{term=}")
    #             ret += term
    #             dingopants = PA.from_expr(term.as_polynomial())
    #             # print(f"{dingopants=}")
    #             dexpr = dexpr - dingopants
    #             dexpr = PA.from_dict({k: v for k, v in dexpr.items() if expand(v) != S.Zero})
    #             # print(f"{dexpr=}")
    #         # print("dingbats")
    #         return ret
    #     except Exception:
    #         raise
    #         #return super().from_expr(expr)

    def mul_expr(self, elem, x):
        if isinstance(x, DomainElement):
            raise TypeError(f"Cannot multiply {type(elem)} with {type(x)}")
        x = sympify(x)
        ind = self.genset.index(x)
        if ind != -1:
            return self.single_variable(elem, ind)
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
        if is_fact_complete_sym(x):
            x = sympify(x)
            if all(a in self.genset for a in genvars(x)) and not any(a in self.genset for a in coeffvars(x)):
                return self.complete_mul(elem, x)
            gens_to_remove = [a for a in genvars(x) if a not in self.genset]

            if any(a in self.genset for a in genvars(x)) and len(gens_to_remove):
                return self.mul_expr(elem, x.split_out_vars(gens_to_remove))

            coeffs_to_remove = [a for a in x.coeff_vars if a in self.genset]

            if len(coeffs_to_remove):
                return self.mul_expr(elem, split_out_vars(x.to_elem_sym(), coeffs_to_remove))

            return self.from_dict({k: (self.handle_sympoly(x)) * v for k, v in elem.items()})
        if isinstance(x, Add):
            return self.sum([self.mul_expr(elem, arg) for arg in x.args])
        if isinstance(x, Mul):
            res = elem
            for arg in x.args:
                res = self.mul_expr(res, arg)
            return res
        if isinstance(x, Pow):
            res = elem
            for _ in range(int(x.args[1])):
                res = self.mul_expr(res, x.args[0])
            return res
        return self.from_dict({k: v * self.domain_new(x) for k, v in elem.items()})

    def new(self, x):
        genset = self.genset
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        if isinstance(x, list) or isinstance(x, tuple):
            p_x = Permutation(x)
            if max([0, *list(p_x.descents())]) > len(self.genset):
                raise NotEnoughGeneratorsError(f"Not enough generators {p_x=} {len(genset)=}")
            elem = self.from_dict({p_x: self.domain.one})
        elif isinstance(x, Permutation):
            if max([0, *list(x.descents())]) > len(self.genset):
                raise NotEnoughGeneratorsError(f"Not enough generators {p_x=} {len(genset)=}")
            elem = self.from_dict({x: self.domain.one})

        elif isinstance(x, DoubleSchubertElement):
            if x.ring.genset == genset:
                return x
            raise ValueError("Different generating set")
        # poly
        # elif isinstance(x, sympy.Poly):
        #     # eject generators we don't want
        #     if not x.has_only_gens():
        #         x, _ = sympy.poly_from_expr(x.as_expr())
        #     new_gens = [g for g in x.gens if self.genset.index(g) != -1]
        #     # end_gens = [g for g in x.gens if self.genset.index(g) == -1]
        #     if len(new_gens) == 0:
        #         # # # logger.debug((f"Didn't find any gens in {x=}")
        #         return self(x.as_expr())
        #     new_gens.sort(key=lambda g: self.genset.index(g))
        #     # expand_gens = [self.genset[i] for i in range(self.genset.index(new_gens[-1])+1)] + end_gens
        #     x = sympy.poly(x, gens=tuple(new_gens))
        #     dct = x.as_dict()
        #     result = 0
        #     for monom, coeff in dct.items():
        #         srt_perm = Permutation.sorting_perm([-i for i in monom])
        #         schub_perm = uncode(sorted(monom, reverse=True))
        #         result += self.from_dict({(schub_perm, 0): coeff}).act(srt_perm)
        #     return result
        else:
            elem = self.from_expr(x)
        return elem


class DoubleSchubertRingDown(DoubleSchubertRing):
    def __hash__(self):
        return hash((self.genset, self.coeff_genset, "fatcabasi"))

    @property
    def double_mul(self):
        return yz.schubmult_double_down

    @property
    def single_mul(self):
        return py.schubmult_py_down

    @cache
    def cached_product(self, u, v, basis2):
        return {k: xreplace_genvars(x, self.coeff_genset, basis2.coeff_genset) for k, x in yz.schubmult_double_down({u: S.One}, v, yz._vars.var_g1, yz._vars.var_g2).items()}

    @cache
    def cached_positive_product(self, u, v, basis2):
        return {k: xreplace_genvars(x, self.coeff_genset, basis2.coeff_genset) for k, x in pos.schubmult_double_down({u: S.One}, v, yz._vars.var_g1, yz._vars.var_g2).items()}

    def printing_term(self, k, prefix="op"):
        return spolymod.DSchubPoly(k, self.genset.label, self.coeff_genset.label, prefix=prefix)


class SingleSchubertRing(DoubleSchubertRing):
    def __init__(self, genset):
        super().__init__(genset, poly_genset(0))
        self.dtype = type("DoubleSchubertElement", (DoubleSchubertElement,), {"ring": self})

    def __str__(self):
        return f"Schubert polynomial ring in {self.genset.label}"

    def __hash__(self):
        return hash((self.genset, self.coeff_genset, "SBS"))

    def _coerce_mul(self, other):
        """Coerce a basis schubert algebra element so it can be multiplied

        Args:
            other (_type_): _description_

        Returns:
            _type_: _description_
        """
        if type(other.ring) is type(self):
            if self.genset == other.ring.genset:
                return other
        if type(other.ring) is DoubleSchubertRing:
            if self.genset == other.ring.genset:
                return other
        if isinstance(other.ring, qsr.QuantumDoubleSchubertRing):
            return self._coerce_mul(other.as_classical())
        return None

    @cache
    def cached_product(self, u, v, basis2):
        if self == basis2:
            return py.schubmult_py({u: S.One}, v)
        return {k: xreplace_genvars(x, poly_genset(0), basis2.coeff_genset if basis2.coeff_genset else poly_genset(0)) for k, x in yz.schubmult_double_pair_generic(u, v).items()}

    @cache
    def cached_positive_product(self, u, v, basis2):
        return self.cached_product(u, v, basis2)

    # need mul sympy

    def single_variable(self, elem, varnum):
        ret = self.zero
        for u, v in elem.items():
            # ret += v * self.domain_new(self.coeff_genset[u[varnum - 1]])
            new_perms = schub_lib.elem_sym_positional_perms(u, 1, varnum)
            for perm, udiff, sign in new_perms:
                # print((perm, udiff, sign))
                if udiff == 1:
                    # print(f"{self.new(perm)=}")
                    ret += self.domain_new(sign * v) * self.new(perm)
                    # print(ret)
        return ret

    def new(self, x):
        genset = self.genset
        if not genset:
            genset = self.genset
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        if isinstance(x, list) or isinstance(x, tuple):
            elem = self.from_dict({Permutation(x): self.domain.one})
        elif isinstance(x, Permutation):
            elem = self.from_dict({x: S.One})
        else:
            elem = self.from_expr(x)
        # print(f"{elem=}")
        # print(f"{x=} {type(x)=} {elem.ring=}")
        return elem

    @property
    def elem_func(self):
        return ElemSym


def DSx(x, genset=GeneratingSet("y"), elem_sym=False, down=False):
    if isinstance(genset, str):
        genset = GeneratingSet(genset)
    if down:
        return DoubleSchubertRingDown(GeneratingSet("x"), genset)(x)
    if elem_sym:
        return ElemDoubleSchubertRing(GeneratingSet("x"), genset)(x)
    return DoubleSchubertRing(GeneratingSet("x"), genset)(x)


Sx = SingleSchubertRing(GeneratingSet("x"))


class ElemDoubleSchubertRing(DoubleSchubertRing):
    def __init__(self, genset, coeff_genset):
        super().__init__(genset, coeff_genset)
        self.dtype = type("DoubleSchubertElement", (DoubleSchubertElement,), {"ring": self})

    def __hash__(self):
        return hash((self.genset, self.coeff_genset, "EDBS"))

    @property
    def replacematch(self):
        def bob(*args, **kwargs):  # noqa: ARG001
            # print(f"{args=} {kwargs=}")
            a = kwargs["a"]
            b = kwargs["b"]
            ind1 = self.genset.index(a)
            ind2 = self.genset.index(b)
            if ind1 != -1:
                if ind2 == -1:
                    return FactorialElemSym(1, 1, [a], [b])
                return FactorialElemSym(1, 1, [a], [self.coeff_genset[1]]) - FactorialElemSym(1, 1, [b], [self.coeff_genset[1]])
            if ind2 != -1:
                return -FactorialElemSym(1, 1, [b], [a])
            if isinstance(a, Symbol) and isinstance(b, Symbol):
                return FactorialElemSym(1, 1, [a], [b])
            return a - b

        return bob

    @property
    def elem_func(self):
        return FactorialElemSym

    def handle_sympoly(self, other):
        return other

    def elem_mul(self, ring_elem, elem):
        indexes = [self.genset.index(a) for a in genvars(elem)]
        ret = self.zero
        elem_sympy = sympify_sympy(elem)
        for k, v in ring_elem.items():
            perm_list = schub_lib.elem_sym_positional_perms(k, degree(elem), *indexes)
            for perm, df, sign in perm_list:
                remaining_vars = [self.coeff_genset[perm[i - 1]] for i in indexes if perm[i - 1] == k[i - 1]]
                coeff = elem_sympy.func(degree(elem) - df, numvars(elem) - df, remaining_vars, coeffvars(elem))
                toadd = self.domain_new(v * sign * coeff) * self(perm)
                # print(f"{toadd=}")
                ret += toadd
        return ret

    def complete_mul(self, elem, x):
        indexes = {self.genset.index(a) for a in genvars(x)}
        ret = self.zero
        x_sympy = sympify_sympy(x)
        for k, v in elem.items():
            perm_list = schub_lib.complete_sym_positional_perms(k, degree(x), *indexes)
            for perm, df, sign in perm_list:
                # print(f"{(perm, df, sign)=}")
                remaining_vars = [self.coeff_genset[perm[i - 1]] for i in {*indexes, *[j + 1 for j in range(len(perm)) if perm[j] != k[j]]}]
                coeff = x_sympy.func(degree(x) - df, numvars(x) + df, remaining_vars, coeffvars(x))  # leave as elem sym
                ret += self.domain_new(sign * v * coeff) * self(perm)
        return ret

    @cache
    def cached_product(self, u, v, basis2):
        return yz.schubmult_double_from_elems({u: self.domain.one}, v, self.coeff_genset, basis2.coeff_genset, elem_func=self.elem_func)

    @cache
    def cached_positive_product(self, u, v, basis2):
        return {k: expand(v) for k, v in yz.schubmult_double_alt_from_elems({u: self.domain.one}, v, self.coeff_genset, basis2.coeff_genset, elem_func=self.elem_func).items()}

    def new(self, x):
        genset = self.genset
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        if isinstance(x, list) or isinstance(x, tuple):
            p_x = Permutation(x)
            if max([0, *list(p_x.descents())]) > len(self.genset):
                raise NotEnoughGeneratorsError(f"Not enough generators {p_x=} {len(genset)=}")
            elem = self.from_dict({p_x: self.domain.one})
        elif isinstance(x, Permutation):
            if max([0, *list(x.descents())]) > len(self.genset):
                raise NotEnoughGeneratorsError(f"Not enough generators {p_x=} {len(genset)=}")
            elem = self.from_dict({x: self.domain.one})
        elif x.ring == self:
            return x
        else:
            elem = self.from_expr(x)
        return elem

    def _coerce_mul(self, other):
        if isinstance(other, BaseSchubertElement):
            if isinstance(other.ring, qsr.QuantumDoubleSchubertRing):
                return other.as_classical()
            if isinstance(other.ring, ElemDoubleSchubertRing):
                return other
            if isinstance(other.ring, DoubleSchubertRing):
                return other
        return None


# class ElemSingleSchubertRing(ElemDoubleSchubertRing):
#     def __init__(self, genset):
#         super().__init__(genset, poly_genset(0))
#         self.dtype = type("DoubleSchubertElement", (DoubleSchubertElement,), {"ring": self})

#     def __hash__(self):
#         return hash((self.genset, self.coeff_genset, "EDBS"))

#     # @property
#     # def replacematch(self):
#     #     def bob(*args, **kwargs):
#     #         # print(f"{args=} {kwargs=}")
#     #         a = kwargs["a"]
#     #         b = kwargs["b"]
#     #         ind1 = self.genset.index(a)
#     #         ind2 = self.genset.index(b)
#     #         if ind1 != -1:
#     #             if ind2 == -1:
#     #                 return FactorialElemSym(1, 1, [a], [b])
#     #             return FactorialElemSym(1, 1, [a], [self.coeff_genset[1]]) - FactorialElemSym(1, 1, [b], [self.coeff_genset[1]])
#     #         if ind2 != -1:
#     #             return -FactorialElemSym(1, 1, [b], [a])
#     #         if isinstance(a, Symbol) and isinstance(b, Symbol):
#     #             return FactorialElemSym(1, 1, [a], [b])
#     #         return a - b

#     #     return bob

#     @property
#     def elem_func(self):
#         return ElemSym

#     def handle_sympoly(self, other):
#         return other

#     def elem_mul(self, ring_elem, elem):
#         indexes = [self.genset.index(a) for a in genvars(elem)]
#         ret = self.zero
#         elem_sympy = sympify_sympy(elem)
#         for k, v in ring_elem.items():
#             perm_list = schub_lib.elem_sym_positional_perms(k, degree(elem), *indexes)
#             for perm, df, sign in perm_list:
#                 remaining_vars = [self.coeff_genset[perm[i - 1]] for i in indexes if perm[i - 1] == k[i - 1]]
#                 coeff = elem_sympy.func(degree(elem) - df, numvars(elem) - df, remaining_vars, coeffvars(elem))
#                 toadd = self.domain_new(v * sign * coeff) * self(perm)
#                 # print(f"{toadd=}")
#                 ret += toadd
#         return ret

#     def complete_mul(self, elem, x):
#         indexes = {self.genset.index(a) for a in genvars(x)}
#         ret = self.zero
#         x_sympy = sympify_sympy(x)
#         for k, v in elem.items():
#             perm_list = schub_lib.complete_sym_positional_perms(k, degree(x), *indexes)
#             for perm, df, sign in perm_list:
#                 # print(f"{(perm, df, sign)=}")
#                 remaining_vars = [self.coeff_genset[perm[i - 1]] for i in {*indexes, *[j + 1 for j in range(len(perm)) if perm[j] != k[j]]}]
#                 coeff = x_sympy.func(degree(x) - df, numvars(x) + df, remaining_vars, coeffvars(x))  # leave as elem sym
#                 ret += self.domain_new(sign * v * coeff) * self(perm)
#         return ret

#     @cache
#     def cached_product(self, u, v, basis2):
#         return yz.schubmult_double_from_elems({u: self.domain.one}, v, self.coeff_genset, basis2.coeff_genset, elem_func=self.elem_func)

#     @cache
#     def cached_positive_product(self, u, v, basis2):
#         return {k: expand(v) for k, v in yz.schubmult_double_alt_from_elems({u: self.domain.one}, v, self.coeff_genset, basis2.coeff_genset, elem_func=self.elem_func).items()}

#     def new(self, x):
#         genset = self.genset
#         if not isinstance(genset, GeneratingSet_base):
#             raise TypeError
#         if isinstance(x, list) or isinstance(x, tuple):
#             p_x = Permutation(x)
#             if max([0, *list(p_x.descents())]) > len(self.genset):
#                 raise NotEnoughGeneratorsError(f"Not enough generators {p_x=} {len(genset)=}")
#             elem = self.from_dict({p_x: self.domain.one})
#         elif isinstance(x, Permutation):
#             if max([0, *list(x.descents())]) > len(self.genset):
#                 raise NotEnoughGeneratorsError(f"Not enough generators {p_x=} {len(genset)=}")
#             elem = self.from_dict({x: self.domain.one})
#         elif x.ring == self:
#             return x
#         else:
#             elem = self.from_expr(x)
#         return elem

#     def _coerce_mul(self, other):
#         if isinstance(other, BaseSchubertElement):
#             if isinstance(other.ring, qsr.QuantumDoubleSchubertRing):
#                 return other.as_classical()
#             if isinstance(other.ring, ElemDoubleSchubertRing):
#                 return other
#             if isinstance(other.ring, DoubleSchubertRing):
#                 return other
#         return None
