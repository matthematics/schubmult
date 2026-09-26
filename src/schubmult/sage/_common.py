"""Shared machinery for Sage parents backed by schubmult rings.

A backed ring is a :class:`~sage.combinat.free_module.CombinatorialFreeModule` indexed by
permutations whose base ring is an infinite polynomial ring over the scalars in the coefficient
alphabets (``y``, ``z``, ``q``, ...). Arithmetic is delegated to a schubmult ring object
(:meth:`SchubmultBackedRing._schub_ring`) and coefficients are converted at the boundary.

Index conventions: Sage variables are 0-indexed (``x0``, ``y_0``, ``q_0``), schubmult's are 1-indexed
(``x_1``, ``y_1``, ``q_1``); see :mod:`schubmult.sage._convert`.
"""

from sage.categories.filtered_algebras_with_basis import FilteredAlgebrasWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.permutation import Permutation, Permutations
from sage.combinat.schubert_polynomial import SchubertPolynomialRing_xbasis
from sage.misc.cachefunc import cached_method
from sage.rings.polynomial.infinite_polynomial_element import InfinitePolynomial
from sage.rings.polynomial.infinite_polynomial_ring import InfinitePolynomialRing
from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ

from ._convert import parse_sage_name, sage_polynomial_to_symengine, symengine_to_sage

X_LETTER = "x"


def to_schubmult_perm(w):
    from schubmult.combinatorics.permutation import Permutation as SPermutation

    return SPermutation(list(w))


def to_sage_perm(w):
    return Permutation(list(w) or [1]).remove_extra_fixed_points()


def genset(letter):
    from schubmult.symbolic.poly.variables import GeneratingSet

    return GeneratingSet(letter)


def coefficient_into(c, T):
    """Move a base-ring coefficient (infinite polynomial in ``a_<i>``) into the finite ring ``T`` with variables ``a<i>``."""
    if not isinstance(c, InfinitePolynomial):
        return T(c)
    p = c.polynomial()
    names = [f"{letter}{idx}" for letter, idx in map(parse_sage_name, p.parent().variable_names())]
    result = T.zero()
    for exps, coeff in p.dict().items():
        term = T(coeff)
        for name, e in zip(names, exps):
            if e:
                term *= T(name) ** e
        result += term
    return result


class SchubmultBackedElement(CombinatorialFreeModule.Element):
    def expand(self):
        r"""
        Expand into a polynomial in ``x0, x1, ...`` and the coefficient variables ``y0, q0, ...``.

        There are `n` variables ``x``, `n` the size of the largest permutation involved (as for
        :meth:`sage.combinat.schubert_polynomial.SchubertPolynomial_class.expand`); coefficient
        letters get exactly the indices that occur.
        """
        import symengine

        P = self.parent()
        R = P.base_ring().base_ring()
        n = max([len(w) for w in self.support()] + [1])
        exprs = {w: symengine.sympify(P._basis_polynomial(w)) for w in self.support()}
        counts = dict.fromkeys(P._alphabets, 0)
        counts[X_LETTER] = n
        for expr in exprs.values():
            for s in expr.free_symbols:
                letter, idx = parse_sage_name(str(s))
                counts[letter] = max(counts.get(letter, 0), idx)  # schubmult index i is Sage index i-1
        for c in self.coefficients():
            if isinstance(c, InfinitePolynomial):
                for v in c.polynomial().variables():
                    letter, idx = parse_sage_name(str(v))
                    counts[letter] = max(counts[letter], idx + 1)
        names = [f"{X_LETTER}{i}" for i in range(counts[X_LETTER])] + [f"{a}{i}" for a in P._alphabets for i in range(counts[a])]
        T = PolynomialRing(R, len(names), names)  # explicit count: one name alone would give a univariate ring
        gens = dict(zip(names, T.gens()))

        def variable(letter, i):
            return gens[f"{letter}{i - 1}"]

        def scalar(q):
            return T(q) if isinstance(q, int) else T(QQ(q.numerator) / QQ(q.denominator))

        result = T.zero()
        for w, c in self:
            result += coefficient_into(c, T) * symengine_to_sage(exprs[w], variable, scalar)
        return result


class SchubmultBackedRing(CombinatorialFreeModule):
    """Base class: subclasses set ``_alphabet`` (basis alphabet letter or ``None``), ``_alphabets`` (all
    coefficient letters, sorted), and implement ``_schub_ring(alphabet)`` and ``_check_basis_perm``.
    """

    Element = SchubmultBackedElement

    def __init__(self, R, alphabets, prefix, name):
        self._alphabets = tuple(alphabets)
        self._name = name
        self._repr_option_bracket = False
        base = InfinitePolynomialRing(R, list(alphabets))
        # filtered, not graded: S_u S_v = sum c^w_{uv} S_w with l(w) <= l(u) + l(v), the coefficients
        # carrying the missing degree (and q terms lower the length further)
        CombinatorialFreeModule.__init__(self, base, Permutations(), category=FilteredAlgebrasWithBasis(base), prefix=prefix)

    def _repr_(self):
        return f"{self._name} over {self.base_ring().base_ring()}"

    # ---- hooks -----------------------------------------------------------------------------

    def _schub_ring(self, alphabet=None):
        """The schubmult ring object computing this basis (with the given second alphabet, if any)."""
        raise NotImplementedError

    def _check_basis_perm(self, w):
        """Raise ``ValueError`` if ``w`` may not index a basis element (e.g. non-parabolic)."""

    def _basis_polynomial(self, w):
        """SymEngine polynomial of the basis element indexed by ``w`` (in schubmult's 1-indexed variables)."""
        return self._schub_ring()(to_schubmult_perm(w)).as_polynomial()

    # ---- conversions -----------------------------------------------------------------------

    def _scalar(self, q):
        B = self.base_ring()
        return B(q) if isinstance(q, int) else B(QQ(q.numerator) / QQ(q.denominator))

    def _variable(self, letter, i):
        return self.base_ring().gen(self._alphabets.index(letter))[i - 1]

    def _convert_dict(self, dct):
        """schubmult ``{Permutation: symengine coeff}`` -> element of ``self``."""
        return self._from_dict({to_sage_perm(w): symengine_to_sage(c, self._variable, self._scalar) for w, c in dct.items()}, remove_zeros=True)

    def _from_polynomial(self, p):
        gensets = {X_LETTER: genset(X_LETTER), **{a: genset(a) for a in self._alphabets}}
        return self._convert_dict(self._schub_ring().from_expr(sage_polynomial_to_symengine(p, gensets)))

    def _from_other_alphabet(self, elem):
        """Expand an element of the same kind of ring with another second alphabet in this basis."""
        ring, other = self._schub_ring(), self._schub_ring(elem.parent()._alphabet)
        one = ring(to_schubmult_perm([]))
        result = self.zero()
        for w, c in elem:
            result += c * self._convert_dict(one * other(to_schubmult_perm(w)))
        return result

    # ---- algebra structure -----------------------------------------------------------------

    @cached_method
    def one_basis(self):
        return self._indices([1])

    def degree_on_basis(self, w):
        return w.length()

    def product_on_basis(self, left, right):
        ring = self._schub_ring()
        return self._convert_dict(ring(to_schubmult_perm(left)) * ring(to_schubmult_perm(right)))

    def _element_constructor_(self, x):
        if isinstance(x, list | tuple):
            x = list(x)
            if x not in Permutations():
                raise ValueError(f"the input {x} is not a valid permutation")
            x = Permutation(x)
        if isinstance(x, Permutation):
            w = x.remove_extra_fixed_points()
            self._check_basis_perm(w)
            return self._from_dict({w: self.base_ring().one()})
        if isinstance(x, Polynomial):  # univariate: re-wrap as a one-variable multivariate polynomial
            x = PolynomialRing(x.base_ring(), 1, x.parent().variable_name())(x)
        if isinstance(x, MPolynomial | InfinitePolynomial):
            return self._from_polynomial(x)
        parent = getattr(x, "parent", None)
        if parent is not None:
            parent = parent()
            if isinstance(parent, SchubertPolynomialRing_xbasis):
                return self._from_polynomial(x.expand())
            if isinstance(parent, SchubmultBackedRing):
                if self._same_kind(parent):
                    if parent._alphabet == self._alphabet:
                        return self._from_dict({w: self.base_ring()(c) for w, c in x})
                    return self._from_other_alphabet(x)
                return self._from_polynomial(x.expand())
        raise TypeError(f"do not know how to make an element of {self} from {x!r}")

    def _same_kind(self, other):
        """Same family of rings (so only the second alphabet may differ)."""
        return type(other) is type(self) and getattr(other, "_parabolic", None) == getattr(self, "_parabolic", None)

    def _coerce_map_from_(self, S):
        if isinstance(S, SchubertPolynomialRing_xbasis):
            return self.base_ring().has_coerce_map_from(S.base_ring())
        if isinstance(S, SchubmultBackedRing):
            if getattr(S, "_parabolic", None) not in (None, getattr(self, "_parabolic", None)):
                return False
            return self.base_ring().has_coerce_map_from(S.base_ring())
        return super()._coerce_map_from_(S)

    def some_elements(self):
        perms = [w for w in ([2, 1], [3, 2, 1], [4, 2, 1, 3], [1, 3, 2], [2, 3, 1], [1, 2, 4, 3]) if self._is_basis_perm(Permutation(w))]
        out = [self.one()]
        if perms:
            out.append(self.one() + 2 * self(perms[0]))
        if len(perms) > 2:
            out.append(self(perms[2]) - self(perms[1]))
        return out

    def _is_basis_perm(self, w):
        try:
            self._check_basis_perm(w)
        except ValueError:
            return False
        return True
