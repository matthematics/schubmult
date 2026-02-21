from functools import cache

from schubmult.combinatorics.permutation import uncode
from schubmult.rings.base_ring import BaseRing, BaseRingElement
from schubmult.symbolic import (
    EXRAW,
    CoercionFailed,
    S,
    expand,
    sstr,
    sympify,
)
from schubmult.utils.logging import get_logger

from ..schubert.schubert_ring import DSx, Sx
from ..schubert.separated_descents import SeparatedDescentsRing
from .free_algebra_basis import SchubertBasis, WordBasis

splugSx = SeparatedDescentsRing(Sx([]).ring)
ADSx = SeparatedDescentsRing(DSx([]).ring)

logger = get_logger(__name__)


# keys are tuples of nonnegative integers
class FreeAlgebraElement(BaseRingElement):

    def poly_inner_product(self, poly, genset, n):
        from schubmult.symbolic.poly.variables import genset_dict_from_expr

        wordish = self.change_basis(WordBasis)
        result = 0

        dct0 = genset_dict_from_expr(poly, genset)
        dct = {}
        for k, v in dct0.items():
            if n is None:
                k = [*k]
                while len(k) > 0 and k[-1] == 0:
                    k.pop()
                dct[tuple(k)] = dct.get(tuple(k), S.Zero) + v
                continue
            if len(k) > n:
                if not all(a == 0 for a in k[n:]):
                    return 0
                dct[k[:n]] = v
            else:
                kpop = (*k, *([0] * (n - len(k))))
                dct[kpop] = v
        # print(f"{dct=}")
        for k, v in wordish.items():
            result += int(v * dct.get(k, S.Zero))
        return result

    def kill_zero(self, fat=False, val=S.Zero):
        spink = self.change_basis(WordBasis)
        spoink = spink.ring.zero
        for k, v in spink.items():
            if fat:
                if 0 in k:
                    v *= val ** k.count(0)
            spoink += v * spink.ring(*[a for a in k if a != 0])
        return spoink.change_basis(self.ring._basis)

    def __hash__(self):
        return hash(set(self.items()))

    def __imul__(self, other):
        return self * other

    def __rmul__(self, other):
        from ..schubert.schubert_ring import DoubleSchubertElement, SingleSchubertRing
        from .free_algebra_basis import SchubertBasis

        if isinstance(other, DoubleSchubertElement):
            if not isinstance(other.ring, SingleSchubertRing):
                other = Sx([]) * other
            ret0 = self.change_basis(SchubertBasis)
            ret = ret0.ring.zero
            for k, v in other.items():
                ret += v * (ret0 / k)
            return ret.change_basis(self.ring._basis)

        try:
            return self.ring.rmul(self, other)
        except CoercionFailed:
            return NotImplemented

    def as_coefficients_dict(self):
        return {self.ring.printing_term(k, self.ring): sympify(v) for k, v in self.items()}

    def expand(self, deep=True, *args, **kwargs):  # noqa: ARG002
        return self.ring.from_dict({k: expand(v, **kwargs) for k, v in self.items()})

    @property
    def free_symbols(self):
        return set()

    def hom_nsym(self):
        # from symengine import Symbol

        fa_elem = self.change_basis(WordBasis)
        # t = Symbol("t")
        t = 1
        ret = fa_elem.ring.zero
        for k, v in fa_elem.items():
            k2 = tuple(a for a in k if a != 0)
            if k2 != k:
                ret += v * (t ** (len(k) - len(k2))) * fa_elem.ring(*k2)
            else:
                ret += v * fa_elem.ring(*k2)
        return ret.change_basis(self.ring._basis)

    @staticmethod
    @cache
    def tup_double_expand(tup):
        res = ADSx([])
        if len(tup) == 0:
            return res
        tup2 = tup[:-1]
        return FreeAlgebraElement.tup_double_expand(tup2) * ADSx(uncode([tup[-1]]), 1)

    @staticmethod
    @cache
    def tup_expand(tup):
        res = splugSx([])
        if len(tup) == 0:
            return res
        if len(tup) == 1:
            return splugSx(uncode(tup), 1)
        mid = len(tup) // 2
        return FreeAlgebraElement.tup_expand(tup[:mid]) * FreeAlgebraElement.tup_expand(tup[mid:])

    def change_basis(self, other_basis):
        new_ring = FreeAlgebra(basis=other_basis)
        ret = new_ring.zero
        tfunc = self.ring._basis.transition(other_basis)
        for k, v in self.items():
            ret += v * new_ring.from_dict(tfunc(k))
        return ret

    def schub_expand(self):
        res = splugSx([]).ring.zero
        for tup, val in self.items():
            res += val * self.__class__.tup_expand(tup)
        return res

    def schub_double_expand(self):
        res = ADSx([]).ring.zero
        for tup, val in self.items():
            res += val * self.__class__.tup_double_expand(tup)
        return res

    def bcoproduct(self):
        T = self.ring @ self.ring
        res = T.zero

        for key, val in self.items():
            res += val * self.ring.bcoproduct_on_basis(key)
        return res

    def factorize(self, j):
        from .free_algebra_basis import FreeAlgebraBasis

        spink = self.change_basis(WordBasis)
        ring = spink.ring @ spink.ring
        ret = ring.zero

        for k, v in spink.items():
            ret += v * ring((k[:j], k[j:]))
        return FreeAlgebraBasis.change_tensor_basis(ret, self.ring._basis, self.ring._basis)

    # def antipode(self):
    #     new_elem = self.change_basis(WordBasis)
    #     ret = new_elem.ring.zero
    #     for k, v in new_elem.items():
    #         to_add = new_elem.ring.one
    #         for a in k:
    #             if a == 0:
    #                 to_add *= (new_elem.ring.one - new_elem.ring(a))
    #             else:
    #                 to_add *= -new_elem.ring(a)
    #         ret += v * to_add
    #         #ret[tuple(reversed(k))] = (S.NegativeOne**(len(k) - k.count(0)))*v
    #     return ret.change_basis(self.ring._basis)

    def __eq__(self, other):
        if isinstance(other, FreeAlgebraElement):
            if other.ring == self.ring:
                diff = self - other
                return all(v == S.Zero for v in diff.values())
        return False
        # return type(self) is type(other) and self.ring == other.ring and dict.__eq__(self, other)

    def __str__(self):
        return sstr(self)

    def _nsymtup(self, tup, R):
        if len(tup) == 0:
            return R.one
        if 0 in tup:
            # return self._nsymtup(tup[:-1], R)
            return R.zero
        return self._nsymtup(tup[:-1], R) * R.from_dict({(tup[-1],): S.One}, R)

    def remove_zeros(self, inserter=S.One):
        new_elem = self.ring.zero
        for k, v in self.items():
            # if 0 in k:
            #     continue
            # add x?
            new_elem += self.ring.from_dict({tuple([a for a in k if a != 0]): v}) * (inserter ** len([a for a in k if a == 0]))
        return new_elem

    def __truediv__(self, other):
        from schubmult.combinatorics.permutation import Permutation

        if isinstance(other, list | tuple | Permutation):
            other = Permutation(other)
            ret = self.ring.zero
            for k, v in self.items():
                ret += v * self.ring.skew_element(k[0], other, k[1])
            return ret
        raise NotImplementedError("Division by non-permutation is not implemented.")

    # def nsymexpand(self):
    #     R = NSym()
    #     ret = R.zero
    #     for k, v in self.items():
    #         ret += v * self._nsymtup(k, R)
    #     return ret

    def split(self, p):
        T = self.ring @ self.ring
        ret = T.zero
        for tup, val in self.items():
            if len(tup) < p:
                ret += val * T((tup, ()))
            else:
                ret += val * T((tup[:p], tup[p:]))
        return ret

    def to_schub(self, sym=False):
        res = Sx([]).ring.zero
        for k, v in self.items():
            res += v * self.ring.tup_to_schub(k, sym=sym)
        return res


class FreeAlgebra(BaseRing):
    def __str__(self):
        return self.__class__.__name__

    def mul_expr(self, elem, x):
        try:
            return self.from_dict({k: x * v for k, v in elem.items()})
        except CoercionFailed:
            # import traceback
            # traceback.print_exc()
            pass
        raise CoercionFailed

    # @cache
    # def j_quasi_recurse(self, alphagod):
    #     from sage.all import ZZ, QuasiSymmetricFunctions

    #     tt = ZZ["t"]
    #     QSym = QuasiSymmetricFunctions(tt)
    #     M = QSym.M()
    #     ret = QSym.zero()
    #     if len(alphagod) == sum(alphagod):
    #         return M[*alphagod]
    #     zero_poly = Sx(uncode(alphagod))
    #     asum = sum(alphagod)
    #     dct = zero_poly.pull_out_gen(zero_poly.ring.genset[1])
    #     for key, _ in dct.items():
    #         oldval = self.j_quasi_recurse(tuple(key.trimcode))
    #         for monom, coeff in dict(oldval).items():
    #             ret += coeff * M[asum - key.inv, *monom]

    #     new_alphagod = [0, *alphagod]
    #     schub_poly = Sx(uncode(new_alphagod))
    #     dct = schub_poly.pull_out_gen(schub_poly.ring.genset[1])
    #     for key, _ in dct.items():
    #         cd = key.trimcode
    #         if len(cd) == 0 or cd[0] == 0 or asum - key.inv == 0 or len(cd) != len(alphagod) or 0 in cd:
    #             continue
    #         oldval = self.j_quasi_recurse(tuple(key.trimcode))
    #         for monom, coeff in dict(oldval).items():
    #             ret += tt.gens()[0] * coeff * M[asum - key.inv, *monom]
    #     return ret

    # @cache
    # def j_quasisymmetric(self, alphagod):
    #     ret0 = self.j_quasi_recurse(alphagod)
    #     from sage.all import ZZ, QuasiSymmetricFunctions

    #     QSym = QuasiSymmetricFunctions(ZZ)
    #     ret = QSym.zero()
    #     M = QSym.M()
    #     for k, v in dict(ret0).items():
    #         ret += int(v.subs(1)) * M[*k]
    #     return ret

    def tensor_schub_expand(self, tensor):
        T = splugSx([]).ring @ splugSx([]).ring
        ret = T.zero
        for (key1, key2), val in tensor.items():
            ret += val * T.ext_multiply(self(key1).schub_expand(), self(key2).schub_expand())
        return ret

    # def tensor_nsym_expand(self, tensor, inserter=S.Zero):
    #     T = NSym() @ NSym()
    #     ret = T.zero
    #     for (key1, key2), val in tensor.items():
    #         ret += val * T.ext_multiply(self(key1).remove_zeros(inserter=inserter).nsymexpand(), self(key2).remove_zeros(inserter=inserter).nsymexpand())
    #     return ret

    def __hash__(self):
        return hash((self.domain, "whatabong", self._basis))

    def __eq__(self, other):
        return type(self) is type(other) and self.domain == other.domain

    def __init__(self, basis=WordBasis, domain=None):
        super().__init__(domain=domain)
        if domain is None:
            self.domain = EXRAW
            self.dom = self.domain
        self._basis = basis
        self.zero_monom = self._basis.zero_monom
        self.dtype = type("FreeAlgebraElement", (FreeAlgebraElement,), {"ring": self})

    @staticmethod
    def right_pad(tup, n):
        if len(tup) < n:
            return (*tup, *([0] * (n - len(tup))))
        return tup

    @cache
    def coproduct_on_basis(self, key):
        T = self @ self
        return T.from_dict(self._basis.coproduct(key))

    @cache
    def bcoproduct_on_basis(self, key):
        T = self @ self
        return T.from_dict(self._basis.bcoproduct(key))

    def _mul_elements(self, elem, other):
        ret = self.zero
        for k0, v0 in elem.items():
            for k, v in other.items():
                ret += self.from_dict(self._basis.product(k0, k, v * v0))
        if self._basis == WordBasis or not FreeAlgebra.CAP:
            return ret
        n = FreeAlgebra.CAP
        ret = {k: v for k, v in ret.items() if len(k[0]) <= n}
        return self.from_dict(ret)

    def from_rc_graph(self, rc_graph):
        return self.from_dict(self._basis.from_rc_graph(rc_graph))

    def matmul(self, elem, other):
        try:
            other = self.domain_new(other)
            return self.from_dict({k: other * v for k, v in elem.items()})
        except Exception:
            pass
        if isinstance(other, FreeAlgebraElement):
            ret = self.zero
            for k0, v0 in elem.items():
                for k, v in other.items():
                    ret += self.from_dict(self._basis.internal_product(k0, k, v * v0))
            return ret
        raise CoercionFailed

    CAP = None

    def new(self, *x):
        if len(x) == 0 and isinstance(x, FreeAlgebraElement):
            return x
        if len(x) == 1 and self._basis.is_key(x[0]):
            return self.from_dict({x[0]: S.One})
        if self._basis.is_key(x):
            return self.from_dict({self._basis.as_key(x): S.One})
        return self.mul_scalar(self.one, x)

    def printing_term(self, k):
        return self._basis.printing_term(k)

    def _coerce_mul(self, other):
        if isinstance(other, FreeAlgebraElement) and other.ring == self:
            return other
        return False

    @property
    def one(self):
        return self.from_dict({self.zero_monom: S.One})

    def from_dict(self, element):
        poly = self.zero
        for monom, coeff in element.items():
            if coeff != S.Zero:
                poly[monom] = coeff
        return poly

    def skew_element(self, w, u, n):
        """Skew schubert by elem sym"""
        return self.from_dict(self._basis.skew_element(w, u, n))

    def domain_new(self, element, orig_domain=None):  # noqa: ARG002
        if isinstance(element, FreeAlgebraElement) or isinstance(element, BaseRingElement):
            raise CoercionFailed("Not a domain element")
        return sympify(element)


FA = FreeAlgebra()

ASx = FreeAlgebra(SchubertBasis)
