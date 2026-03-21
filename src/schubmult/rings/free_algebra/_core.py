from functools import cache
from itertools import zip_longest

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
    """Element of a free algebra, stored as a dict mapping basis keys to coefficients.

    Keys are tuples of nonnegative integers (words in the word basis) or
    basis-specific keys depending on the parent ring's basis. Supports
    arithmetic operations, basis changes, and word-level operations like
    injection, prefix, suffix, and interval extraction.
    """

    def interleave(self, other, zero_pad=True):
        """Interleave two elements letter-by-letter in the word basis.

        Converts both elements to WordBasis, then interleaves each pair of
        words by alternating entries (a1, b1, a2, b2, ...). Shorter words
        are zero-padded when ``zero_pad`` is True.

        Args:
            other: Another FreeAlgebraElement to interleave with.
            zero_pad: If True, pad shorter words with zeros.

        Returns:
            The interleaved element in the original basis.
        """
        if not isinstance(other, FreeAlgebraElement):
            return NotImplemented

        self_word = self.change_basis(WordBasis)
        other_word = other.change_basis(WordBasis)

        r = self_word.ring.zero
        for word, coeff in self_word.items():
            for word2, coeff2 in other_word.items():
                interleaved = []
                for a, b in zip_longest(word, word2, fillvalue=None):
                    #if a is not None:
                    if zero_pad or len(a) == len(b):
                        interleaved.append(a if a is not None else 0)
                        interleaved.append(b if b is not None else 0)
                r += coeff * coeff2 * self_word.ring(*interleaved)

        return r.change_basis(self.ring._basis)

    def inject(self, i, other):
        """Insert another element's words at position *i* in this element's words.

        Delegates to the current basis's ``inject`` classmethod.

        Args:
            i: Nonnegative integer insertion index.
            other: Another FreeAlgebraElement to inject.

        Returns:
            A new element with the injected words.
        """
        if not isinstance(other, FreeAlgebraElement):
            return NotImplemented
        if not isinstance(i, int):
            raise TypeError(f"Expected integer index, got {type(i)}")
        if i < 0:
            raise IndexError(f"Insertion index {i} must be nonnegative")

        result = self.ring.zero
        for k0, v0 in self.items():
            for k, v in other.items():
                result += self.ring.from_dict(self.ring._basis.inject(k0, i, k, v * v0))
        return result

    def prefix(self, length):
        """Extract the first *length* letters of each word.

        Delegates to the current basis's ``prefix`` classmethod.

        Args:
            length: Nonnegative integer prefix length.

        Returns:
            A new element containing the prefixes.
        """
        if not isinstance(length, int):
            raise TypeError(f"Expected integer length, got {type(length)}")
        if length < 0:
            raise IndexError(f"Prefix length {length} must be nonnegative")

        result = self.ring.zero
        for k, v in self.items():
            result += self.ring.from_dict(self.ring._basis.prefix(k, length, v))
        return result

    def suffix(self, length):
        """Extract the last *length* letters of each word.

        Delegates to the current basis's ``suffix`` classmethod.

        Args:
            length: Nonnegative integer suffix length.

        Returns:
            A new element containing the suffixes.
        """
        if not isinstance(length, int):
            raise TypeError(f"Expected integer length, got {type(length)}")
        if length < 0:
            raise IndexError(f"Suffix length {length} must be nonnegative")

        result = self.ring.zero
        for k, v in self.items():
            result += self.ring.from_dict(self.ring._basis.suffix(k, length, v))
        return result

    def interval(self, start, stop):
        """Extract a subword from position *start* to *stop* in each word.

        Delegates to the current basis's ``interval`` classmethod.

        Args:
            start: Nonnegative start index (inclusive).
            stop: Stop index (exclusive).

        Returns:
            A new element containing the subwords.
        """
        if not isinstance(start, int) or not isinstance(stop, int):
            raise TypeError(f"Expected integer indices, got {type(start)} and {type(stop)}")
        if start < 0:
            raise IndexError(f"Start index {start} must be nonnegative")

        result = self.ring.zero
        for k, v in self.items():
            result += self.ring.from_dict(self.ring._basis.interval(k, start, stop, v))
        return result

    def poly_inner_product(self, poly, genset, n):
        """Compute the inner product of this element with a polynomial.

        Converts to WordBasis and pairs coefficient-by-coefficient with
        the monomial expansion of *poly* in *genset*.

        Args:
            poly: A polynomial expression.
            genset: The generating set of variables for the polynomial.
            n: Number of variables to use (or None for automatic).

        Returns:
            The integer inner product value.
        """
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
        """Remove zeros from each word key.

        In the word basis, strips all zero entries from each key. When
        *fat* is True, multiplies the coefficient by ``val`` raised to
        the number of removed zeros instead of simply dropping them.

        Args:
            fat: If True, weight by ``val`` per removed zero.
            val: The value to raise per zero when *fat* is True.

        Returns:
            A new element with zeros removed, in the original basis.
        """
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
        """Return a dict mapping printing terms to sympified coefficients."""
        return {self.ring.printing_term(k, self.ring): sympify(v) for k, v in self.items()}

    def expand(self, deep=True, *args, **kwargs):  # noqa: ARG002
        """Expand all coefficients symbolically."""
        return self.ring.from_dict({k: expand(v, **kwargs) for k, v in self.items()})

    def pairing(self, other):
        """Compute the pairing of this element with *other* via the monomial basis.

        Converts *self* to WordBasis and *other* to MonomialBasis, then
        sums products of matching coefficients.

        Args:
            other: Another element to pair with.

        Returns:
            The integer pairing value.
        """
        from schubmult.rings.polynomial_algebra import MonomialBasis
        monomo = self.change_basis(WordBasis)
        monomo2 = other.change_basis(MonomialBasis)
        result = 0
        for k, v in monomo.items():
            result += int(v * monomo2.get(k, S.Zero))
        return result

    @property
    def free_symbols(self):
        return set()

    def hom_nsym(self):
        """Apply the homomorphism to noncommutative symmetric functions.

        Strips zero entries from each word key (dropping them from the word)
        and returns the result in the original basis.
        """
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
        """Expand a word tuple into the double Schubert basis via Pieri products."""
        res = ADSx([])
        if len(tup) == 0:
            return res
        tup2 = tup[:-1]
        return FreeAlgebraElement.tup_double_expand(tup2) * ADSx(uncode([tup[-1]]), 1)

    @staticmethod
    @cache
    def tup_expand(tup):
        """Expand a word tuple into the single Schubert basis via divide-and-conquer Pieri products."""
        res = splugSx([])
        if len(tup) == 0:
            return res
        if len(tup) == 1:
            return splugSx(uncode(tup), 1)
        mid = len(tup) // 2
        return FreeAlgebraElement.tup_expand(tup[:mid]) * FreeAlgebraElement.tup_expand(tup[mid:])

    def change_basis(self, other_basis):
        """Convert this element to another basis.

        Args:
            other_basis: The target basis class (e.g. WordBasis, SchubertBasis).

        Returns:
            A new FreeAlgebraElement in the target basis's ring.
        """
        new_ring = FreeAlgebra(basis=other_basis)
        ret = new_ring.zero
        tfunc = self.ring._basis.transition(other_basis)
        for k, v in self.items():
            ret += v * new_ring.from_dict(tfunc(k))
        return ret

    def schub_expand(self):
        """Expand this element into a single Schubert polynomial ring element."""
        res = splugSx([]).ring.zero
        for tup, val in self.items():
            res += val * self.__class__.tup_expand(tup)
        return res

    def schub_double_expand(self):
        """Expand this element into a double Schubert polynomial ring element."""
        res = ADSx([]).ring.zero
        for tup, val in self.items():
            res += val * self.__class__.tup_double_expand(tup)
        return res

    def bcoproduct(self):
        """Compute the bar-coproduct of this element in the tensor ring."""
        T = self.ring @ self.ring
        res = T.zero

        for key, val in self.items():
            res += val * self.ring.bcoproduct_on_basis(key)
        return res

    def factorize(self, j):
        """Split each word at position *j*, returning a tensor element.

        In the word basis, each word ``w`` maps to ``(w[:j], w[j:])``.
        The result is expressed in the tensor ring of the original basis.

        Args:
            j: Position at which to split each word.

        Returns:
            An element of the tensor product ring.
        """
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
        """Remove zero entries from each key, weighting by *inserter* per zero removed.

        Args:
            inserter: Scalar multiplied per removed zero (default 1).

        Returns:
            A new element with zeros stripped from keys.
        """
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
        """Split each word at position *p* into a tensor element.

        Words shorter than *p* are placed entirely in the left factor.

        Args:
            p: Position at which to split.

        Returns:
            An element of the tensor product ring.
        """
        T = self.ring @ self.ring
        ret = T.zero
        for tup, val in self.items():
            if len(tup) < p:
                ret += val * T((tup, ()))
            else:
                ret += val * T((tup[:p], tup[p:]))
        return ret

    def to_schub(self, sym=False):
        """Convert to a Schubert ring element via ``tup_to_schub``.

        Args:
            sym: If True, use symmetric expansion.

        Returns:
            An element of the single Schubert ring.
        """
        res = Sx([]).ring.zero
        for k, v in self.items():
            res += v * self.ring.tup_to_schub(k, sym=sym)
        return res


class FreeAlgebra(BaseRing):
    """Free algebra ring with a configurable basis.

    The algebra operates on :class:`FreeAlgebraElement` instances whose keys
    are determined by the chosen basis (default :class:`WordBasis`). Supports
    multiplication, tensor products, coproducts, and basis changes.

    Args:
        basis: The basis class to use (default ``WordBasis``).
        domain: Coefficient domain (default ``EXRAW``).
    """

    def __str__(self):
        return self.__class__.__name__

    def mul_expr(self, elem, x):
        """Multiply every coefficient of *elem* by the scalar *x*."""
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
        """Expand a tensor element into the Schubert polynomial tensor ring."""
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
        """Initialize a FreeAlgebra with the given basis and coefficient domain."""
        super().__init__(domain=domain)
        if domain is None:
            self.domain = EXRAW
            self.dom = self.domain
        self._basis = basis
        self.zero_monom = self._basis.zero_monom
        self.dtype = type("FreeAlgebraElement", (FreeAlgebraElement,), {"ring": self})

    @staticmethod
    def right_pad(tup, n):
        """Right-pad *tup* with zeros to length *n*."""
        if len(tup) < n:
            return (*tup, *([0] * (n - len(tup))))
        return tup

    @cache
    def coproduct_on_basis(self, key):
        """Compute the coproduct of a single basis key in the tensor ring."""
        T = self @ self
        return T.from_dict(self._basis.coproduct(key))

    @cache
    def bcoproduct_on_basis(self, key):
        """Compute the bar-coproduct of a single basis key in the tensor ring."""
        T = self @ self
        return T.from_dict(self._basis.bcoproduct(key))

    def mul(self, elem, other):
        """Multiply two elements via the basis product rule."""
        try:
            ret = self.zero
            for k0, v0 in elem.items():
                for k, v in other.items():
                    ret += self.from_dict(self._basis.product(k0, k, v * v0))
            if self._basis == WordBasis or not FreeAlgebra.CAP:
                return ret
            n = FreeAlgebra.CAP
            ret = {k: v for k, v in ret.items() if len(k[0]) <= n}
            return self.from_dict(ret)
        except Exception:
            return super().mul(elem, other)

    def from_rc_graph(self, rc_graph):
        """Create an element from an RC graph."""
        return self.from_dict(self._basis.from_rc_graph(rc_graph))

    def matmul(self, elem, other):
        """Internal product (``@`` operator) or scalar multiplication.

        If *other* is a scalar, multiplies all coefficients. If *other* is
        a FreeAlgebraElement, computes the internal product via the basis.
        """
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
        """Create a new element from the given key or arguments."""
        if len(x) == 0 and isinstance(x, FreeAlgebraElement):
            return x
        if len(x) == 1 and self._basis.is_key(x[0]):
            return self.from_dict({x[0]: S.One})
        if self._basis.is_key(x):
            return self.from_dict({self._basis.as_key(x): S.One})
        return self.mul_scalar(self.one, x)

    def printing_term(self, k):
        """Return the display symbol for basis key *k*."""
        return self._basis.printing_term(k)

    def _coerce_mul(self, other):
        if isinstance(other, FreeAlgebraElement) and other.ring == self:
            return other
        return False

    @property
    def one(self):
        return self.from_dict({self.zero_monom: S.One})

    def from_dict(self, element):
        """Construct an element from a dict of ``{key: coefficient}`` pairs."""
        poly = self.zero
        for monom, coeff in element.items():
            if coeff != S.Zero:
                poly[monom] = coeff
        return poly

    def skew_element(self, w, u, n):
        """Skew schubert by elem sym"""
        return self.from_dict(self._basis.skew_element(w, u, n))

    def domain_new(self, element, orig_domain=None):  # noqa: ARG002
        """Coerce a raw value into the coefficient domain."""
        if isinstance(element, FreeAlgebraElement) or isinstance(element, BaseRingElement):
            raise CoercionFailed("Not a domain element")
        return sympify(element)


FA = FreeAlgebra()

ASx = FreeAlgebra(SchubertBasis)
