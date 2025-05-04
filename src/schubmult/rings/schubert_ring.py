from functools import cache, cached_property

import symengine
import sympy
from symengine import Add, Mul, Pow, S
from sympy import CoercionFailed

import schubmult.rings.quantum_schubert_ring as qsr
import schubmult.schub_lib.double as yz
import schubmult.schub_lib.schub_lib as schub_lib
import schubmult.schub_lib.single as py
import schubmult.utils.ring_utils as utils
from schubmult.perm_lib import Permutation, uncode
from schubmult.utils.logging import get_logger
from schubmult.utils.perm_utils import add_perm_dict

from .abstract_schub_poly import AbstractSchubPoly
from .backend import CompleteSym, ElemSym, sympify
from .base_schubert_ring import BaseSchubertElement, BaseSchubertRing
from .poly_lib import elem_sym_poly, xreplace_genvars
from .schub_poly import schubpoly_classical_from_elems
from .symmetric_polynomials.sympy.functions import split_out_vars
from .tensor_ring import TensorRing
from .variables import CustomGeneratingSet, GeneratingSet, GeneratingSet_base, MaskedGeneratingSet

_pretty_schub_char = "𝔖"  # noqa: RUF001

logger = get_logger(__name__)


class DoubleSchubertElement(BaseSchubertElement):
    """Algebra with sympy coefficients
    and a dict basis
    """

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
                    result += self.ring.from_dict({new_perm: v}).mult_poly(sympy.prod([(new - coeff_gens[index2]) for index2 in index_list]))
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
                    for index2 in index_list:
                        toadd *= self.ring.genset[index2] - gen
                    ret += toadd * val * new_basis(~new_perm)
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

    # def elem_coproduct(self, *indices):
    #     genset = self.ring.genset
    #     genset = self.ring.coeff_genset
    #     gname1 = f"{genset.label}_A"  # "("+", ".join([f"{genset.label}_{i}" for i in indices])+")"
    #     gname2 = f"{genset.label}_B"  # f"{genset.label}\\{{"+", ".join([f"{genset.label}_{i}" for i in indices])+"}"
    #     gens2 = MaskedGeneratingSet(genset, indices)
    #     gens1 = gens2.complement()
    #     gens1.set_label(gname1)
    #     gens2.set_label(gname2)
    #     R = TensorRing(DoubleSchubertRing(gens1),DoubleSchubertRing(gens2))
    #     def elem_sym(p, k, varl1, varl2):

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


class DSchubPoly(AbstractSchubPoly):
    is_Atom = True

    def __new__(cls, k, basis):
        return DSchubPoly.__xnew_cached__(cls, k, basis)

    @staticmethod
    def __xnew__(_class, k, basis):
        return AbstractSchubPoly.__new__(_class, k, basis)

    def _sympystr(self, printer):
        key = self._key
        if self._key == Permutation([]):
            return printer.doprint(1)
        if self.ring.coeff_genset.label is None:
            return printer.doprint(f"S{self.ring.genset.label}({printer.doprint(key)})")
        return printer.doprint(f"DS{self.ring.genset.label}({printer.doprint(key)}, {self.ring.coeff_genset.label})")

    def _pretty(self, printer):
        key = self._key
        gl = self.ring.genset.label
        if key == Permutation([]):
            return printer._print(1)
        subscript = printer._print(int("".join([str(i) for i in key])))
        if self.ring.coeff_genset.label is None:
            return printer._print_Function(sympy.Function(f"{_pretty_schub_char}_{subscript}")(sympy.Symbol(gl)))
        return printer._print_Function(sympy.Function(f"{_pretty_schub_char}_{subscript}")(sympy.Symbol(f"{self.ring.genset.label}; {self.ring.coeff_genset.label}")))

    def _latex(self, printer):
        key = self._key
        gl = self.ring.genset.label
        if key == Permutation([]):
            return printer._print(1)
        subscript = printer._print(key)
        if self.ring.coeff_genset.label is None:
            return printer._print_Function(sympy.Function("\\mathfrak{S}" + f"_{'{' + subscript + '}'}")(sympy.Symbol(gl)))
        return printer._print_Function(sympy.Function("\\mathfrak{S}" + f"_{'{' + subscript + '}'}")(sympy.Symbol(f"{{{self.ring.genset.label}}}; {{{self.ring.coeff_genset.label}}}")))

    @staticmethod
    @cache
    def __xnew_cached__(_class, k, basis):
        return DSchubPoly.__xnew__(_class, k, basis)


class DoubleSchubertRing(BaseSchubertRing):
    def __hash__(self):
        return hash((self.genset, self.coeff_genset, "DBS"))

    def __init__(self, genset, coeff_genset, domain=None):
        super().__init__(genset, coeff_genset, domain)
        self.dtype = type("DoubleSchubertElement", (DoubleSchubertElement,), {"ring": self})

    def __str__(self):
        return f"Double Schubert polynomial ring in {self.genset.label} and {self.coeff_genset.label}"

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

    #     return self.mul_sympy(elem, other)

    def printing_term(self, k):
        return DSchubPoly(k, self)

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
        genset = self.genset

        class esf(sympy.Function):
            @classmethod
            def eval(cls, *x):
                pass

            def _eval_expand_func(self, **_):
                if len(self.args) == 2:
                    return elem_sym_poly(int(self.args[0]), int(self.args[1]), genset[1:], utils.poly_ring(0))
                return elem_sym_poly(int(self.args[0]), int(self.args[1]), genset[1:], self.args[2:])

        return esf

    @property
    def elem_mul_type(self):
        return ElemSym

    # specifically symengine
    def elem_mul(self, ring_elem, elem):
        indexes = [self.genset.index(a) for a in elem.genvars]
        ret = self.zero
        for k, v in ring_elem.items():
            perm_list = schub_lib.elem_sym_positional_perms(k, elem._p, *indexes)
            for perm, df, sign in perm_list:
                remaining_vars = [self.coeff_genset[perm[i - 1]] for i in indexes if perm[i - 1] == k[i - 1]]
                coeff = ElemSym(elem._p - df, elem._k - df, remaining_vars, elem.coeffvars)  # leave as elem sym
                ret += self.domain_new(v * sympy.expand_func(coeff)) * self(perm)
        return ret

    @property
    def symbol_elem_func(self):
        return ElemSym

    def elem_sym_subs(self, kk):
        elems = []
        for k in range(1, kk + 1):
            for p in range(1, k + 1):
                elems += [(sympy.Symbol(f"e_{p}_{k}"), elem_sym_poly(p, k, self.genset[1:], utils.poly_ring(0)))]
        return dict(elems)

    @staticmethod
    def flip(elem):
        R = DoubleSchubertRing(CustomGeneratingSet([0, *elem.coeffvars]), CustomGeneratingSet([0, *elem.genvars]))
        p = elem._p
        K = elem._k + 1 - p
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
        return {k: xreplace_genvars(x, self.coeff_genset, basis2.coeff_genset if basis2.coeff_genset else utils.poly_ring(0)) for k, x in yz.schubmult_double_pair_generic(u, v).items()}

    @cache
    def cached_positive_product(self, u, v, basis2):
        return {k: xreplace_genvars(x, self.coeff_genset, basis2.coeff_genset if basis2.coeff_genset else utils.poly_ring(0)) for k, x in yz.schubmult_generic_partial_posify(u, v).items()}

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
        return self.from_dict({(schub_perm, utils.NoneVar): S.One}).act(srt_perm)

    @cache
    def cached_schubpoly(self, k):
        return schubpoly_classical_from_elems(k, self.genset, self.coeff_genset, elem_func=elem_sym_poly)

    def complete_mul(self, elem, x):
        indexes = {self.genset.index(a) for a in x.genvars}
        ret = self.zero
        for k, v in elem.items():
            perm_list = schub_lib.complete_sym_positional_perms(k, x._p, *indexes)
            for perm, df, sign in perm_list:
                remaining_vars = [self.coeff_genset[perm[i - 1]] for i in {*indexes, *[j + 1 for j in range(len(perm)) if perm[j] != k[j]]}]
                coeff = CompleteSym(x._p - df, x._k + df, remaining_vars, x.coeffvars)  # leave as elem sym
                ret += self.domain_new(sign * v * sympy.expand_func(coeff)) * self(perm)
        return ret

    def handle_sympoly(self, other):
        return sympy.expand_func(other)

    def single_variable(self, elem, varnum):
        ret = self.zero
        for u, v in elem.items():
            ret += v * self.domain_new(self.coeff_genset[u[varnum - 1]]) * self(u)
            new_perms = schub_lib.elem_sym_positional_perms(u, 1, varnum)
            for perm, udiff, sign in new_perms:
                if udiff == 1:
                    ret += self.domain_new(sign * v) * self(perm)
        return ret

    def mul_sympy(self, elem, x):
        _Add = Add
        _Mul = Mul
        _Pow = Pow
        try:
            x = symengine.sympify(x)
        except Exception:
            x = sympy.sympify(x)
            _Add = sympy.Add
            _Mul = sympy.Mul
            _Pow = sympy.Pow
        ind = self.genset.index(x)
        if ind != -1:
            return self.single_variable(elem, ind)
        if isinstance(sympify(x), ElemSym):
            if all(self.genset.index(a)!=-1 for a in x.genvars) and not any(self.genset.index(a)!=-1 for a in x.coeffvars):
                return self.elem_mul(elem, x)

            gens_to_remove = [a for a in x.genvars if a not in self.gensWet]
            if any(self.genset.index(a)!=-1 for a in x.genvars) and len(gens_to_remove):
                return self.mul_sympy(elem, x.split_out_vars(gens_to_remove))

            coeffs_to_remove = [a for a in x.coeffvars if a in self.genset]

            if any(a in self.genset for a in x.coeffvars) and len(coeffs_to_remove):
                return self.mul_sympy(elem.split_out_vars(x.to_complete_sym(), coeffs_to_remove))
            return self.from_dict({k: self.domain_new(self.handle_sympoly(x)) * v for k, v in elem.items()})
        if isinstance(sympify(x), CompleteSym):
            if all(a in self.genset for a in x.genvars) and not any(a in self.genset for a in x.coeffvars):
                return self.complete_mul(elem, x)
            gens_to_remove = [a for a in x.genvars if a not in self.genset]

            if any(a in self.genset for a in x.genvars) and len(gens_to_remove):
                return self.mul_sympy(elem, x.split_out_vars(gens_to_remove))

            coeffs_to_remove = [a for a in x.coeff_vars if a in self.genset]

            if len(coeffs_to_remove):
                return self.mul_sympy(elem, split_out_vars(x.to_elem_sym(), coeffs_to_remove))

            return self.from_dict({k: self.domain_new(self.handle_sympoly(x)) * v for k, v in elem.items()})
        if isinstance(x, _Add):
            return self.sum([self.mul_sympy(elem, arg) for arg in x.args])
        if isinstance(x, _Mul):
            res = elem
            for arg in x.args:
                res = self.mul_sympy(res, arg)
            return res
        if isinstance(x, _Pow):
            res = elem
            for _ in range(int(x.args[1])):
                res = self.mul_sympy(res, x.args[0])
            return res
        return self.from_dict({k: v * self.domain_new(x) for k, v in elem.items()})

    
    def new(self, x):
        genset = self.genset
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        if isinstance(x, list) or isinstance(x, tuple):
            p_x = Permutation(x)
            if max([0, *list(p_x.descents())]) > len(self.genset):
                raise utils.NotEnoughGeneratorsError(f"Not enough generators {p_x=} {len(genset)=}")
            elem = self.from_dict({p_x: self.domain.one})
        elif isinstance(x, Permutation):
            if max([0, *list(x.descents())]) > len(self.genset):
                raise utils.NotEnoughGeneratorsError(f"Not enough generators {p_x=} {len(genset)=}")
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
        #         result += self.from_dict({(schub_perm, utils.NoneVar): coeff}).act(srt_perm)
        #     return result
        else:
            elem = self.from_sympy(x)
        return elem


def DSx(x, genset=GeneratingSet("y"), elem_sym=False):
    if isinstance(genset, str):
        genset = GeneratingSet(genset)
    if elem_sym:
        return ElemDoubleSchubertRing(GeneratingSet("x"), genset)(x)
    return DoubleSchubertRing(GeneratingSet("x"), genset)(x)


class SingleSchubertRing(DoubleSchubertRing):
    def __init__(self, genset):
        super().__init__(genset, utils.poly_ring(utils.NoneVar))

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
        return {k: xreplace_genvars(x, utils.poly_ring(0), basis2.coeff_genset if basis2.coeff_genset else utils.poly_ring(0)) for k, x in yz.schubmult_double_pair_generic(u, v).items()}

    @cache
    def cached_positive_product(self, u, v, basis2):
        return self.cached_product(u, v, basis2)

    # need mul sympy

    # def single_variable(self, elem, varnum):
    #     ret = self.zero
    #     for u, v in elem.items():
    #         # ret += v * self.domain_new(self.coeff_genset[u[varnum - 1]])
    #         new_perms = schub_lib.elem_sym_positional_perms(u, 1, varnum)
    #         for perm, udiff, sign in new_perms:
    #             if udiff == 1:
    #                 ret += self.domain_new(sign * v) * self(perm)
    #     return ret

    def mul_sympy(self, elem, x):
        _Add = Add
        _Mul = Mul
        _Pow = Pow
        try:
            x = symengine.sympify(x)
        except Exception:
            x = sympy.sympify(x)
            _Add = sympy.Add
            _Mul = sympy.Mul
            _Pow = sympy.Pow
        ind = self.genset.index(x)
        if ind != -1:
            return self.single_variable(elem, ind)
        if isinstance(sympify(x), ElemSym):
            if all(a in self.genset for a in x.genvars) and not any(a in self.genset for a in x.coeffvars):
                return self.elem_mul(elem, x)

            gens_to_remove = [a for a in x.genvars if a not in self.gensWet]
            if any(a in self.genset for a in x.genvars) and len(gens_to_remove):
                return self.mul_sympy(elem, x.split_out_vars(gens_to_remove))

            coeffs_to_remove = [a for a in x.coeffvars if a in self.genset]

            if any(a in self.genset for a in x.coeffvars) and len(coeffs_to_remove):
                return self.mul_sympy(elem.split_out_vars(x.to_complete_sym(), coeffs_to_remove))
            return self.from_dict({k: self.domain_new(self.handle_sympoly(x)) * v for k, v in elem.items()})
        if isinstance(sympify(x), CompleteSym):
            if all(a in self.genset for a in x.genvars) and not any(a in self.genset for a in x.coeffvars):
                return self.complete_mul(elem, x)
            gens_to_remove = [a for a in x.genvars if a not in self.genset]

            if any(a in self.genset for a in x.genvars) and len(gens_to_remove):
                return self.mul_sympy(elem, x.split_out_vars(gens_to_remove))

            coeffs_to_remove = [a for a in x.coeff_vars if a in self.genset]

            if len(coeffs_to_remove):
                return self.mul_sympy(elem, split_out_vars(x.to_elem_sym(), coeffs_to_remove))

            return self.from_dict({k: self.domain_new(self.handle_sympoly(x)) * v for k, v in elem.items()})
        if isinstance(x, _Add):
            return self.sum([self.mul_sympy(elem, arg) for arg in x.args])
        if isinstance(x, _Mul):
            res = elem
            for arg in x.args:
                res = self.mul_sympy(res, arg)
            return res
        if isinstance(x, _Pow):
            res = elem
            for _ in range(int(x.args[1])):
                res = self.mul_sympy(res, x.args[0])
            return res
        return self.from_dict({k: v * self.domain_new(x) for k, v in elem.items()})

    # def from_sympy(self, x):
    #     if isinstance(x, BaseSchubertElement):
    #         if x.ring == self:
    #             return x
    #     x = sympify(x)
    #     ind = self.genset.index(x)
    #     if ind != -1:
    #         return self.from_dict(py.mult_poly_py({Permutation([]): self.domain.one}, x, self.genset))
    #     if isinstance(x, Add):
    #         return self.sum([self.from_sympy(arg) for arg in x.args])
    #     if isinstance(x, Mul):
    #         res = self.one
    #         for arg in x.args:
    #             res *= self.from_sympy(arg)
    #         return res
    #     if isinstance(x, Pow):
    #         return self.from_sympy(x.args[0]) ** int(x.args[1])
    #     return self.from_dict({Permutation([]): self.domain_new(x)})

    # def from_sympy(self, x):
    #     x = sympify(x)
    #     result = py.mult_poly_py({Permutation([]): 1}, x, self.genset)
    #     return self.from_dict(result)

    def new(self, x):
        genset = self.genset
        if not genset:
            genset = self.genset
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        if isinstance(x, list) or isinstance(x, tuple):
            elem = self.from_dict({Permutation(x): self.domain.one})
        elif isinstance(x, Permutation):
            elem = self.from_dict({x: 1})
        elif isinstance(x, DoubleSchubertElement):
            if x.ring.genset == genset:
                elem = DoubleSchubertElement(x, self)  # , self)
            else:
                return self(x.expand())
        else:
            elem = self.from_sympy(x)
        return elem


Sx = SingleSchubertRing(GeneratingSet("x"))

# class ElemSymRingElem(Expr):
#     def

# class ExprWrapped(sympy.Basic):

#     # @property
#     # def __slots__(self):
#     #     return (*self._old_slots, "domain")

#     @property
#     def domain(self):
#         return self.args[0]

#     @property
#     def args(self):
#         return self._args[1:]

#     def __new__(cls, *args, domain=None, _class=sympy.Basic):
#         obj = _.__new__(cls, *args)
#         return obj

# class ElemSymDomain(EXRAW):


#     def new(self, a):
#         a = sympify(a)
#         return a


# ESDOM = ElemSymDomain()


class ElemDoubleSchubertRing(DoubleSchubertRing):
    def __init__(self, genset, coeff_genset):
        super().__init__(genset, coeff_genset)
        self.dtype = type("DoubleSchubertElement", (DoubleSchubertElement,), {"ring": self})

    _elem_sym_pattern = sympy.Wild("a") - sympy.Wild("b")
    _elem_sym_pattern2 = sympy.Wild("a") + sympy.Wild("b")

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
                    return ElemSym(1, 1, [a], [b])
                return ElemSym(1, 1, [a], [self.coeff_genset[1]]) - ElemSym(1, 1, [b], [self.coeff_genset[1]])
            if ind2 != -1:
                return -ElemSym(1, 1, [b], [a])
            if isinstance(a, sympy.Symbol) and isinstance(b, sympy.Symbol):
                return ElemSym(1, 1, [a], [b])
            return a - b

        return bob

    # def domain_new(self, element, orig_domain=None):
    #     elem = self.domain.new(element)
    #     if element.has_free(*self.symbols):
    #         raise CoercionFailed(f"{element} contains an element of the set of generators")
    #     return elem

    # def domain_new(self, element, orig_domain=None):  # noqa: ARG002
    #     element = sympify(element)
    #     if not element.has_free(*self.symbols):
    #         return element
    #     raise CoercionFailed(f"{element} contains an element of the set of generators")

    @property
    def elem_func(self):
        return ElemSym

    def handle_sympoly(self, other):
        return other

    # def single_variable(self, elem, varnum):
    #     ret = self.zero
    #     for u, v in elem:
    #         ret += v * self.domain_new(ElemSym(1,1,[self.coeff_genset[u[varnum-1]])
    #         new_perms  = schub_lib.elem_sym_positional_perms(u, 1, varnum)
    #         for perm, udiff, sign in new_perms:
    #             if udiff == 1:
    #                 ret += self.domain_new(sign * v) * self(perm)
    #     return ret

    def elem_mul(self, ring_elem, elem):
        if not all(a in self.genset for a in elem.genvars):
            gens_to_remove = [a for a in elem.genvars if a not in self.genset]
            return ring_elem * elem.split_out_vars(gens_to_remove)
        indexes = [self.genset.index(a) for a in elem.genvars]
        ret = self.zero
        for k, v in ring_elem.items():
            perm_list = schub_lib.elem_sym_positional_perms(k, elem._p, *indexes)
            for perm, df, sign in perm_list:
                remaining_vars = [self.coeff_genset[perm[i - 1]] for i in indexes if perm[i - 1] == k[i - 1]]
                coeff = ElemSym(elem._p - df, elem._k - df, remaining_vars, elem.coeffvars)
                toadd = self.domain_new(v * sign * coeff) * self(perm)
                # print(f"{toadd=}")
                ret += toadd
        return ret

    def complete_mul(self, elem, x):
        indexes = {self.genset.index(a) for a in x.genvars}
        ret = self.zero
        for k, v in elem.items():
            perm_list = schub_lib.complete_sym_positional_perms(k, x._p, *indexes)
            for perm, df, sign in perm_list:
                # print(f"{(perm, df, sign)=}")
                remaining_vars = [self.coeff_genset[perm[i - 1]] for i in {*indexes, *[j + 1 for j in range(len(perm)) if perm[j] != k[j]]}]
                coeff = CompleteSym(x._p - df, x._k + df, remaining_vars, x.coeffvars)  # leave as elem sym
                ret += self.domain_new(sign * v * coeff) * self(perm)
        return ret

    @cache
    def cached_product(self, u, v, basis2):
        return yz.schubmult_double_from_elems({u: 1}, v, self.coeff_genset, basis2.coeff_genset, elem_func=self.elem_func)

    @cache
    def cached_positive_product(self, u, v, basis2):
        return self.cached_product(u, v, basis2)

    def new(self, x):
        genset = self.genset
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        if isinstance(x, list) or isinstance(x, tuple):
            p_x = Permutation(x)
            if max([0, *list(p_x.descents())]) > len(self.genset):
                raise utils.NotEnoughGeneratorsError(f"Not enough generators {p_x=} {len(genset)=}")
            elem = self.from_dict({p_x: self.domain.one})
        elif isinstance(x, Permutation):
            if max([0, *list(x.descents())]) > len(self.genset):
                raise utils.NotEnoughGeneratorsError(f"Not enough generators {p_x=} {len(genset)=}")
            elem = self.from_dict({x: self.domain.one})
        elif isinstance(x, DoubleSchubertElement):
            if x.ring.genset == genset:
                return x
            raise ValueError("Different generating set")
        else:
            elem = self.from_sympy(x)
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
