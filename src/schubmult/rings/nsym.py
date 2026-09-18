"""`NSym`: a `FreeAlgebra` indexed by compositions, multiplied through the separated-descents Schubert ring.

A key is a composition ``alpha`` (a tuple of positive integers), printed ``N(alpha)``, and is
identified with the Schubert key ``(uncode(alpha - 1), len(alpha))`` of
`schubmult.rings.schubert.separated_descents.SeparatedDescentsRing` via `NSym.sepify` /
`NSym.from_sep`. The product is the separated-descents product transported back to
compositions; when ``FreeAlgebra.CAP`` is set the result is truncated to keys of at most that
length. Right multiplication by a Schubert element acts by the skew operation ``/``.
"""

from schubmult.combinatorics.permutation import uncode
from schubmult.symbolic import (
    EXRAW,
    Add,
    CoercionFailed,
    S,
    sstr,
    sympify,
    sympify_sympy,
    sympy_Add,
    sympy_Mul,
)

from .free_algebra import FreeAlgebra, FreeAlgebraElement
from .printing import GenericPrintingTerm
from .schubert.schubert_ring import Sx
from .schubert.separated_descents import SeparatedDescentsRing

splugSx = SeparatedDescentsRing(Sx([]).ring)


class NSym(FreeAlgebra):
    """Free algebra on compositions with the separated-descents product. See the module docstring."""

    def __hash__(self):
        return hash((self.domain, "whatabong2"))

    def __init__(self, domain=None):
        """Create the ring over ``domain`` (default ``EXRAW``); the empty composition is the identity."""
        if domain:
            self.domain = domain
        else:
            self.domain = EXRAW
        self.dom = self.domain
        self.zero_monom = ()
        self.dtype = type("NSymElement", (NSymElement,), {"ring": self})

    def printing_term(self, k):
        """Display as ``N(alpha)``."""
        return GenericPrintingTerm(k, "N")

    def rmul(self, elem, other):
        """Scale coefficients by the scalar ``other``."""
        # print(f"{self=} {elem=} {other=}")
        if isinstance(other, NSymElement):
            raise NotImplementedError
        return self.from_dict({k: v * other for k, v in elem.items()})

    def sepify(self, elem):
        """Map ``alpha -> (uncode(alpha - 1), len(alpha))`` into the separated-descents Schubert ring."""
        return splugSx([]).ring.from_dict({(uncode([a - 1 for a in k]), len(k)): v for k, v in elem.items()})

    def from_sep(self, elem):
        """Inverse of `sepify`: pad or cut the code of ``perm`` to length ``n`` and add 1 to each entry."""
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
        """Scalar multiplication, or the separated-descents product of two elements (truncated by
        ``FreeAlgebra.CAP`` if set).
        """
        # print(f"{self=} {elem=} {other=}")
        try:
            other = self.domain_new(other)
            # print(f"{other=} {type(other)=}")
            return self.from_dict({k: other * v for k, v in elem.items()})
        except CoercionFailed:
            pass
        if isinstance(other, NSymElement):
            if not FreeAlgebra.CAP:
                return self.from_sep(self.sepify(elem) * self.sepify(other))
            n = FreeAlgebra.CAP
            dd = self.from_sep(self.sepify(elem) * self.sepify(other))
            return {k: v for k, v in dd.items() if len(k[0]) <= n}
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
    """Element of `NSym`: a dict from compositions to coefficients with SymPy-compatible printing."""

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
        """Scalar on the left, or a Schubert element acting by the skew operation ``self / perm``."""
        from .schubert.schubert_ring import DoubleSchubertElement, SingleSchubertRing

        if isinstance(other, DoubleSchubertElement):
            if not isinstance(other.ring, SingleSchubertRing):
                other = Sx([]) * other
            # ret0 = self.change_basis(Schubert)
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
