from functools import cache

from schubmult.perm_lib import uncode
from schubmult.symbolic import (
    EXRAW,
    Add,
    CoercionFailed,
    CompositeDomain,
    DefaultPrinting,
    DomainElement,
    Ring,
    S,
    expand,
    sstr,
    sympify,
    sympify_sympy,
    sympy_Add,
    sympy_Mul,
)
from schubmult.utils.logging import get_logger
from schubmult.utils.perm_utils import add_perm_dict

from .abstract_schub_poly import GenericPrintingTerm
from .base_schubert_ring import BaseSchubertElement
from .free_algebra_basis import WordBasis
from .schubert_ring import DSx, Sx
from .separated_descents import SeparatedDescentsRing

splugSx = SeparatedDescentsRing(Sx([]).ring)
ADSx = SeparatedDescentsRing(DSx([]).ring)

logger = get_logger(__name__)


# keys are tuples of nonnegative integers
class FreeAlgebraElement(DomainElement, DefaultPrinting, dict):
    precedence = 40

    __sympy__ = True

    def parent(self):
        return self.ring

    def poly_inner_product(self, poly, genset, n):
        from schubmult.rings.variables import genset_dict_from_expr
        from schubmult.symbolic import prod, sympify

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

    def eval(self, *args):
        pass

    def __hash__(self):
        return hash(set(self.items()))

    def _sympystr(self, printer):
        if len(self.keys()) == 0:
            return printer._print(S.Zero)
        if printer.order in ("old", "none"):  # needed to avoid infinite recursion
            return printer._print_Add(sympy_Add(*self.as_ordered_terms()), order="lex")
        return printer._print_Add(sympy_Add(*self.as_ordered_terms()))

    def _pretty(self, printer):
        if len(self.keys()) == 0:
            return printer._print(S.Zero)
        if printer.order in ("old", "none"):  # needed to avoid infinite recursion
            return printer._print_Add(self, order="lex")
        return printer._print_Add(sympy_Add(*self.as_ordered_terms()))

    def _latex(self, printer):
        if len(self.keys()) == 0:
            return printer._print(S.Zero)
        if printer.order in ("old", "none"):  # needed to avoid infinite recursion
            return printer._print_Add(self, order="lex")
        return printer._print_Add(sympy_Add(*self.as_ordered_terms()))

    def as_terms(self):
        if len(self.keys()) == 0:
            return [sympify_sympy(S.Zero)]
        return [self[k] if k == () else sympy_Mul(sympify_sympy(self[k]), self.ring.printing_term(k)) for k in self.keys()]

    def as_ordered_terms(self, *_, **__):
        if len(self.keys()) == 0:
            return [sympify(S.Zero)]
        return [((self[k]) if k == () else sympy_Mul(sympify_sympy(self[k]), self.ring.printing_term(k))) for k in sorted(self.keys())]

    def __add__(self, other):
        if isinstance(other, FreeAlgebraElement):
            if self.ring == other.ring:
                return self.ring.add(self, other)
            return other.__radd__(self)
        try:
            other = self.ring.domain_new(other)
            other = self.ring.from_dict({(): other})
            return self.ring.add(self, other)
        except CoercionFailed:
            pass
        try:
            new_other = self.ring(other)
            return self.__add__(new_other)
        except CoercionFailed:
            return other.__radd__(self)

    def __radd__(self, other):
        try:
            other = self.ring.domain_new(other)
            other = self.ring.from_dict({(): other})
            return self.ring.add(other, self)
        except CoercionFailed:
            pass
        try:
            new_other = self.ring(other)
            return new_other.__add__(self)
        except CoercionFailed:
            return NotImplemented

    def __sub__(self, other):
        if isinstance(other, FreeAlgebraElement):
            if self.ring == other.ring:
                return self.ring.sub(self, other)
            return other.__rsub__(self)
        try:
            other = self.ring.domain_new(other)
            other = self.ring.from_dict({(): other})
            return self.ring.sub(self, other)
        except CoercionFailed:
            pass
        try:
            new_other = self.ring(other)
            return self.__sub__(new_other)
        except CoercionFailed:
            return other.__rsub__(self)

    def __rsub__(self, other):
        try:
            other = self.ring.domain_new(other)
            other = self.ring.from_dict({(): other})
            return self.ring.sub(other, self)
        except CoercionFailed:
            pass
        try:
            new_other = self.ring(other)
            return new_other.__sub__(self)
        except CoercionFailed:
            return NotImplemented

    def __neg__(self):
        return self.ring.neg(self)

    def __mul__(self, other):
        try:
            return self.ring.mul(self, other)
        except CoercionFailed:
            return other.__rmul__(self)

    def __matmul__(self, other):
        return self.ring.matmul(self, other)

    def __rmul__(self, other):
        from .schubert_ring import DoubleSchubertElement, SingleSchubertRing

        if isinstance(other, DoubleSchubertElement):
            if not isinstance(other.ring, SingleSchubertRing):
                other = Sx([]) * other
            ret = self.ring.zero
            for k, v in other.items():
                ret += v * (self / k)
            return ret

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

    def coproduct(self):
        T = self.ring @ self.ring
        res = T.zero

        for key, val in self.items():
            res += val * self.ring.coproduct_on_basis(key)
        return res

    def bcoproduct(self):
        T = self.ring @ self.ring
        res = T.zero

        for key, val in self.items():
            res += val * self.ring.bcoproduct_on_basis(key)
        return res

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

    def as_expr(self):
        return Add(*self.as_terms())

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
        from schubmult.perm_lib import Permutation

        if isinstance(other, list | tuple | Permutation):
            other = Permutation(other)
            ret = self.ring.zero
            for k, v in self.items():
                ret += v * self.ring.skew_element(k[0], other, k[1])
            return ret
        raise NotImplementedError("Division by non-permutation is not implemented.")

    def nsymexpand(self):
        R = NSym()
        ret = R.zero
        for k, v in self.items():
            ret += v * self._nsymtup(k, R)
        return ret

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


class FreeAlgebra(Ring, CompositeDomain):
    def __str__(self):
        return self.__class__.__name__

    def tensor_schub_expand(self, tensor):
        T = splugSx([]).ring @ splugSx([]).ring
        ret = T.zero
        for (key1, key2), val in tensor.items():
            ret += val * T.ext_multiply(self(key1).schub_expand(), self(key2).schub_expand())
        return ret

    def tensor_nsym_expand(self, tensor, inserter=S.Zero):
        T = NSym() @ NSym()
        ret = T.zero
        for (key1, key2), val in tensor.items():
            ret += val * T.ext_multiply(self(key1).remove_zeros(inserter=inserter).nsymexpand(), self(key2).remove_zeros(inserter=inserter).nsymexpand())
        return ret

    def __hash__(self):
        return hash((self.domain, "whatabong", self._basis))

    def __eq__(self, other):
        return type(self) is type(other) and self.domain == other.domain

    def __init__(self, basis=WordBasis, domain=None):
        if domain:
            self.domain = domain
        else:
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

    def __matmul__(self, other):
        from .tensor_ring import TensorRing

        return TensorRing(self, other)

    @cache
    def coproduct_on_basis(self, key):
        T = self @ self
        return T.from_dict(self._basis.coproduct(key))

    @cache
    def bcoproduct_on_basis(self, key):
        T = self @ self
        return T.from_dict(self._basis.bcoproduct(key))

    def add(self, elem, other):
        return self.from_dict(add_perm_dict(elem, other))

    def sub(self, elem, other):
        return self.from_dict(add_perm_dict(elem, {k: -v for k, v in other.items()}))

    def neg(self, elem):
        return self.from_dict({k: -v for k, v in elem.items()})

    def rmul(self, elem, other):
        # print(f"{self=} {elem=} {other=}")
        if isinstance(other, FreeAlgebraElement):
            raise NotImplementedError
        return self.from_dict({k: v * other for k, v in elem.items()})

    def mul(self, elem, other):
        try:
            other = self.domain_new(other)
            return self.from_dict({k: other * v for k, v in elem.items()})
        except Exception:
            pass
        if isinstance(other, FreeAlgebraElement):
            ret = self.zero
            for k0, v0 in elem.items():
                for k, v in other.items():
                    ret += self.from_dict(self._basis.product(k0, k, v * v0))
            return ret
        raise CoercionFailed

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

    def to_domain(self):
        return self

    def new(self, *x):
        if len(x) == 0 and isinstance(x, FreeAlgebraElement):
            return x
        if self._basis.is_key(x):
            return self.from_dict({self._basis.as_key(x): S.One})
        return self.mul_scalar(self.one, x)

    def printing_term(self, k):
        return self._basis.printing_term(k)

    def _coerce_mul(self, other): ...

    @property
    def one(self):
        return self.from_dict({(): S.One})

    def from_dict(self, element):
        poly = self.zero
        for monom, coeff in element.items():
            if coeff != self.domain.zero:
                poly[monom] = coeff
        return poly

    def skew_element(self, w, u, n):
        """Skew schubert by elem sym"""
        return self.from_dict(self._basis.skew_element(w, u, n))

    @property
    def zero(self):
        return self.dtype()

    def domain_new(self, element, orig_domain=None):  # noqa: ARG002
        if isinstance(element, FreeAlgebraElement) or isinstance(element, BaseSchubertElement):
            raise CoercionFailed("Not a domain element")
        return sympify(element)


FA = FreeAlgebra()


class NSym(FreeAlgebra):
    def __hash__(self):
        return hash((self.domain, "whatabong2"))

    def __init__(self, domain=None):
        if domain:
            self.domain = domain
        else:
            self.domain = EXRAW
        self.dom = self.domain
        self.zero_monom = ()
        self.dtype = type("NSymElement", (NSymElement,), {"ring": self})

    def printing_term(self, k):
        return GenericPrintingTerm(k, "N")

    def rmul(self, elem, other):
        # print(f"{self=} {elem=} {other=}")
        if isinstance(other, NSymElement):
            raise NotImplementedError
        return self.from_dict({k: v * other for k, v in elem.items()})

    def sepify(self, elem):
        return splugSx([]).ring.from_dict({(uncode([a - 1 for a in k]), len(k)): v for k, v in elem.items()})

    def from_sep(self, elem):
        dct = {}
        for (k, n), v in elem.items():
            cd = k.code
            if len(cd) < n:
                cd += [0] * (n - len(cd))
            elif len(cd) > n:
                cd = cd[:n]
            cd = tuple([c + 1 for c in cd])
            dct[cd] = v
        return self.from_dict(dct)

    def mul(self, elem, other):
        # print(f"{self=} {elem=} {other=}")
        try:
            other = self.domain_new(other)
            # print(f"{other=} {type(other)=}")
            return self.from_dict({k: other * v for k, v in elem.items()})
        except Exception:
            pass
        if isinstance(other, NSymElement):
            return self.from_sep(self.sepify(elem) * self.sepify(other))
            # ret = self.zero
            # for k0, v0 in elem.items():
            #     for k, v in other.items():
            #         if len(k0) > 0 and len(k) > 0:
            #             new_key0 = (*k0[:-1], k0[-1] + k[0], *k[1:])
            #             ret += self.from_dict({new_key0: v * v0})
            #         new_key1 = (*k0, *k)
            #         ret += self.from_dict({new_key1: v * v0})
            # return ret
        raise CoercionFailed


class NSymElement(FreeAlgebraElement):
    precedence = 40

    __sympy__ = True

    def parent(self):
        return self.ring

    def eval(self, *args):
        pass

    def _sympystr(self, printer):
        if len(self.keys()) == 0:
            return printer._print(S.Zero)
        if printer.order in ("old", "none"):  # needed to avoid infinite recursion
            return printer._print_Add(sympy_Add(*self.as_ordered_terms()), order="lex")
        return printer._print_Add(sympy_Add(*self.as_ordered_terms()))

    def _pretty(self, printer):
        if len(self.keys()) == 0:
            return printer._print(S.Zero)
        if printer.order in ("old", "none"):  # needed to avoid infinite recursion
            return printer._print_Add(self, order="lex")
        return printer._print_Add(sympy_Add(*self.as_ordered_terms()))

    def _latex(self, printer):
        if len(self.keys()) == 0:
            return printer._print(S.Zero)
        if printer.order in ("old", "none"):  # needed to avoid infinite recursion
            return printer._print_Add(self, order="lex")
        return printer._print_Add(sympy_Add(*self.as_ordered_terms()))

    def as_terms(self):
        if len(self.keys()) == 0:
            return [sympify_sympy(S.Zero)]
        return [self[k] if k == () else sympy_Mul(sympify_sympy(self[k]), self.ring.printing_term(k)) for k in self.keys()]

    def as_ordered_terms(self, *_, **__):
        if len(self.keys()) == 0:
            return [sympify(S.Zero)]
        return [((self[k]) if k == () else sympy_Mul(sympify_sympy(self[k]), self.ring.printing_term(k))) for k in sorted(self.keys())]

    def __add__(self, other):
        if isinstance(other, NSymElement):
            if self.ring == other.ring:
                return self.ring.add(self, other)
            return other.__radd__(self)
        try:
            other = self.ring.domain_new(other)
            other = self.ring.from_dict({(): other})
            return self.ring.add(self, other)
        except CoercionFailed:
            pass
        try:
            new_other = self.ring(other)
            return self.__add__(new_other)
        except CoercionFailed:
            return other.__radd__(self)

    def __radd__(self, other):
        try:
            other = self.ring.domain_new(other)
            other = self.ring.from_dict({(): other})
            return self.ring.add(other, self)
        except CoercionFailed:
            pass
        try:
            new_other = self.ring(other)
            return new_other.__add__(self)
        except CoercionFailed:
            return NotImplemented

    def __sub__(self, other):
        if isinstance(other, NSymElement):
            if self.ring == other.ring:
                return self.ring.sub(self, other)
            return other.__rsub__(self)
        try:
            other = self.ring.domain_new(other)
            other = self.ring.from_dict({(): other})
            return self.ring.sub(self, other)
        except CoercionFailed:
            pass
        try:
            new_other = self.ring(other)
            return self.__sub__(new_other)
        except CoercionFailed:
            return other.__rsub__(self)

    def __rsub__(self, other):
        try:
            other = self.ring.domain_new(other)
            other = self.ring.from_dict({(): other})
            return self.ring.sub(other, self)
        except CoercionFailed:
            pass
        try:
            new_other = self.ring(other)
            return new_other.__sub__(self)
        except CoercionFailed:
            return NotImplemented

    def __neg__(self):
        return self.ring.neg(self)

    def __mul__(self, other):
        try:
            return self.ring.mul(self, other)
        except CoercionFailed:
            return other.__rmul__(self)

    def __rmul__(self, other):
        from .schubert_ring import DoubleSchubertElement, SingleSchubertRing

        if isinstance(other, DoubleSchubertElement):
            if not isinstance(other.ring, SingleSchubertRing):
                other = Sx([]) * other
            ret = self.ring.zero
            for k, v in other.items():
                ret += v * (self / k)
            return ret

        try:
            return self.ring.rmul(self, other)
        except CoercionFailed:
            return NotImplemented

    def as_coefficients_dict(self):
        return {self.ring.printing_term(k, self.ring): sympify(v) for k, v in self.items()}

    @property
    def free_symbols(self):
        return set()

    def as_expr(self):
        return Add(*self.as_terms())

    def __eq__(self, other):
        return type(self) is type(other) and self.ring == other.ring and dict.__eq__(self, other)

    def __str__(self):
        return sstr(self)
