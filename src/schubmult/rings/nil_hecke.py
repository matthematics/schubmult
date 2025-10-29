from sympy import Expr

from schubmult.mult.double import schubmult_double_down
from schubmult.schub_lib.perm_lib import Permutation
from schubmult.symbolic import (
    EXRAW,
    Add,
    CoercionFailed,
    CompositeDomain,
    DefaultPrinting,
    DomainElement,
    Mul,
    Pow,
    Ring,
    S,
    Symbol,
    expand,
    prod,
    sstr,
    sympify,
    sympify_sympy,
    sympy_Add,
    sympy_Mul,
)
from schubmult.utils.logging import get_logger
from schubmult.utils.perm_utils import add_perm_dict

from .base_schubert_ring import BaseSchubertElement
from .schubert_ring import DoubleSchubertElement, SingleSchubertRing
from .variables import GeneratingSet

logger = get_logger(__name__)


class NilHeckeElement(DomainElement, DefaultPrinting, dict):
    _op_priority = 1e200
    precedence = 40

    __sympy__ = True

    def apply(self, other):
        if not isinstance(other, DoubleSchubertElement):
            raise NotImplementedError
        ret = other.ring.zero
        for k0, v0 in self.items():
            for k, v in other.items():
                perm = k * (~k0)
                if perm.inv == k.inv - k0.inv:
                    new_elem = other.ring.from_dict({perm: v})
                    ret += v0 * new_elem
        return ret

    def __reduce__(self):
        return (self.__class__, self.items())

    def parent(self):
        return self.ring

    def has_free(self, *args):
        return any(s in args for s in self.free_symbols)

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
        return [((self[k]) if k == Permutation([]) else sympy_Mul(sympify_sympy(self[k]), self.ring.printing_term(k))) for k in self.keys()]

    def as_ordered_terms(self, *_, **__):
        if len(self.keys()) == 0:
            return [sympify(S.Zero)]
        return [((self[k]) if k == Permutation([]) else sympy_Mul(sympify_sympy(self[k]), self.ring.printing_term(k))) for k in sorted(self.keys(), key=lambda kk: (kk.inv, tuple(kk)))]

    def __add__(self, other):
        if isinstance(other, NilHeckeElement):
            if self.ring == other.ring:
                return self.ring.add(self, other)
            return other.__radd__(self)
        try:
            other = self.ring.domain_new(other)
            other = self.ring.from_dict({Permutation([]): other})
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
            other = self.ring.from_dict({Permutation([]): other})
            return self.ring.add(other, self)
        except CoercionFailed:
            pass
        try:
            new_other = self.ring(other)
            return new_other.__add__(self)
        except CoercionFailed:
            return NotImplemented

    def __sub__(self, other):
        if isinstance(other, NilHeckeElement):
            if self.ring == other.ring:
                return self.ring.sub(self, other)
            return other.__rsub__(self)
        try:
            other = self.ring.domain_new(other)
            other = self.ring.from_dict({Permutation([]): other})
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
            other = self.ring.from_dict({Permutation([]): other})
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
        try:
            return self.ring.rmul(self, other)
        except CoercionFailed:
            return NotImplemented

    def as_coefficients_dict(self):
        return {self.ring.printing_term(k, self.ring): sympify(v) for k, v in self.items()}

    def expand(self, deep=True, *args, **kwargs):  # noqa: ARG002
        return self.ring.from_dict({k: expand(v, **kwargs) for k, v in self.items()})

    def as_expr(self):
        return Add(*self.as_terms())

    def __eq__(self, other):
        return type(self) is type(other) and self.ring == other.ring and dict.__eq__(self, other)

    def __str__(self):
        return sstr(self)


class NilHeckeRing(Ring, CompositeDomain):
    def __str__(self):
        return self.__class__.__name__

    def __eq__(self, other):
        return type(self) is type(other) and self.genset == other.genset

    def to_sympy(self, elem):
        return elem.as_expr()

    def isobaric(self, perm):
        perm = Permutation(perm)
        if perm.inv == 0:
            return self.one
        index = max(perm.descents())
        elem = self.isobaric(perm.swap(index, index + 1))
        elem = elem * self.from_dict({Permutation([]).swap(index, index + 1): S.One})
        return self.mul(elem, self.genset[index + 1])

    def fgp_operator(self, k, length, q_var=GeneratingSet("q")):
        # xk
        elem = self.zero
        elem += self.genset[k] * self.one
        for i in range(1, k):
            perm = Permutation([])
            perm = perm.swap(i - 1, k - 1)
            elem += self.from_dict({perm: -prod([q_var[j] for j in range(i, k)])})
        for i in range(k + 1, length + 1):
            perm = Permutation([])
            perm = perm.swap(k - 1, i - 1)
            elem += self.from_dict({perm: prod([q_var[j] for j in range(k, i)])})
        return elem

    def subs_fgp(self, poly, length):
        if self.genset.index(poly) != -1:
            return self.fgp_operator(self.genset.index(poly), length)
        if isinstance(poly, Add):
            return self.sum([self.subs_fgp(arg, length) for arg in poly.args])
        if isinstance(poly, Mul):
            return prod([self.subs_fgp(arg, length) for arg in poly.args])
        if isinstance(poly, Pow):
            return prod([self.subs_fgp(poly.args[0], length)] * int(poly.args[1]))
        return poly

    def __init__(self, genset, domain=None):
        self._genset = genset
        self._sring = SingleSchubertRing(genset)
        self.symbols = list(genset)
        if domain:
            self.domain = domain
        else:
            self.domain = EXRAW
        self.dom = self.domain
        self.zero_monom = Permutation([])
        self.dtype = type("NilHeckeElement", (NilHeckeElement,), {"ring": self})

    def add(self, elem, other):
        # print(f"{elem=}")
        # print(f"{other=}")
        return self.from_dict(add_perm_dict(elem, other))

    def sub(self, elem, other):
        return self.from_dict(add_perm_dict(elem, {k: -v for k, v in other.items()}))

    def neg(self, elem):
        return self.from_dict({k: -v for k, v in elem.items()})

    def mul_scalar(self, elem, other):
        try:
            other = self.domain_new(other)
            # print(f"{other=} {type(other)=}")
            return self.from_dict({k: other * v for k, v in elem.items()})
        except Exception:
            pass
        mul_elem = other
        if isinstance(other, BaseSchubertElement):
            if not isinstance(other, DoubleSchubertElement) or not other.ring.genset == self.genset:
                mul_elem = self._sring.from_expr(other.as_polynomial())
        else:
            mul_elem = self._sring.from_expr(other)

        ret = self.zero
        for k, v in mul_elem.items():
            ret += self.from_dict({k2: v * v2 for k2, v2 in schubmult_double_down(elem, k, self.genset, mul_elem.ring.coeff_genset).items()})
        return ret

    def mul_perm(self, elem, perm):
        dct = {}
        for k, v in elem.items():
            newperm = k * perm
            if newperm.inv == k.inv + perm.inv:
                dct[newperm] = v
        return self.from_dict(dct)

    def rmul(self, elem, other):
        # print(f"{self=} {elem=} {other=}")
        if isinstance(other, NilHeckeElement):
            raise NotImplementedError
        if isinstance(other, BaseSchubertElement):
            other = other.as_polynomial()
        return self.from_dict({k: v * other for k, v in elem.items()})

    def mul(self, elem, other):
        # print("Fry the leg")
        # print(f"{self=} {elem=} {other=}")
        try:
            other = self.domain_new(other)
            # print(f"{other=} {type(other)=}")
            return self.from_dict({k: other * v for k, v in elem.items()})
        except Exception:
            pass
        # print("FasOjfao")
        if isinstance(other, NilHeckeElement):
            ret = self.zero
            for k, v in other.items():
                new_elem = self.mul_perm(self.mul_scalar(elem, v), k)
                ret += new_elem
            return ret
        return self.mul_scalar(elem, other)

    def to_domain(self):
        return self

    def from_sympy(self, expr):
        return self.from_expr(expr)

    def new(self, x):
        if isinstance(x, NilHeckeElement):
            return x
        if isinstance(x, list) or isinstance(x, tuple):
            return self.from_dict({Permutation(x): S.One})
        if isinstance(x, Permutation):
            return self.from_dict({x: S.One})
        return self.mul_scalar(self.one, x)

    def __hash__(self):
        return hash(self.genset)

    def printing_term(self, k):
        # from sympy import Symbol

        # return Symbol(f"df({k})", commutative=False)
        class NilHeckeTerm(Expr):
            def _sympystr(self, printer):
                return printer._print(f"df({sstr(k)})")

            def _pretty(self, printer):
                return printer._print(f"\u2202({sstr(k)})")

            def _latex(self, printer):
                return printer._print_Symbol(Symbol(f"\\partial^{sstr(k)}"))

        return NilHeckeTerm()

    def _coerce_mul(self, other): ...

    @property
    def one(self):
        return self.from_dict({Permutation([]): S.One})

    def from_dict(self, element, orig_domain=None):  # noqa: ARG002
        # print(f"{element=}")

        poly = self.zero

        for monom, coeff in element.items():
            if coeff != self.domain.zero:
                poly[monom] = coeff
        return poly

    @property
    def zero(self):
        return self.dtype()

    def domain_new(self, element, orig_domain=None):  # noqa: ARG002
        # print(f"{element=} {type(element)=} bagels {type(sympify(element))=} {sympify(element).has(*self.symbols)=}")
        if isinstance(element, NilHeckeElement) or isinstance(element, BaseSchubertElement):
            raise CoercionFailed("Not a domain element")
        if not any(x in self.genset for x in sympify_sympy(element).free_symbols):
            return sympify(element)
        raise CoercionFailed(f"{element} contains an element of the set of generators")

    @property
    def genset(self):
        return self._genset

    def from_expr(self, x):
        return self.mul_scalar(self.one, x)


df = NilHeckeRing(GeneratingSet("x"))
