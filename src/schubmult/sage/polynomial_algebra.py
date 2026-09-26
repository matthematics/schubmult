r"""
The polynomial ring with its combinatorial bases

:func:`PolynomialAlgebra` is the polynomial ring `R[x_0, x_1, \ldots]` (or `R[x_0, \ldots, x_{n-1}]`)
as a parent with realizations (in the sense of
:class:`sage.categories.with_realizations.WithRealizations`, like :class:`SymmetricFunctions`): one
:class:`CombinatorialFreeModule` per basis, with coercions between them. The bases and their
conversions are computed by schubmult's :class:`~schubmult.rings.polynomial_algebra.PolynomialAlgebra`
bases and Schubert kernels.

Each basis is indexed the way its polynomials are indexed in the literature:

- ``schubert``, ``grothendieck``: by permutations, `\mathfrak S_w` and the `\beta = 1` Grothendieck
  polynomials `\mathfrak G_w` (prefixes ``S``, ``G``);
- ``monomial``, ``key``, ``fundamental_slide``, ``monomial_slide``, ``forest``, ``glide``, ``lascoux``,
  ``grove``: by weak compositions `\alpha`: `x^\alpha`, the key polynomials `\kappa_\alpha`, the slide
  polynomials of Assaf-Searles, the forest polynomials of Nadeau-Spink-Tewari, and the K-theoretic
  (`\beta = 1`) glide, Lascoux and grove polynomials (prefixes ``x``, ``k``, ``F``, ``M``, ``P``, ``Gl``,
  ``L``, ``Gr``);
- ``elementary``: products `\prod_{j < n} e_{a_j}(x_0, \ldots, x_{j-1}) \cdot \prod_i e_{b_i}(x_0, \ldots, x_{n-1})`
  of elementary symmetric polynomials, a basis of the polynomial ring in `n` variables *for each fixed
  `n`* -- so it is only available in ``PolynomialAlgebra(R, n)`` (prefix ``E``; see :meth:`PolynomialAlgebra.elementary`).

schubmult's ``PolynomialAlgebra`` is graded by the number of variables (a key of length `n` lives in
the `n`-variable slice and keys of different lengths multiply to zero -- the structure dual to the
free algebra). ``PolynomialAlgebra(R)`` is the plain polynomial ring instead: trailing zeros of a
composition do not matter, `\mathfrak S_w` is `\mathfrak S_w` however many variables are in play, and
products are polynomial products; the layer pads keys to a common number of variables before calling
schubmult. ``PolynomialAlgebra(R, n)`` is the `n`-variable slice itself, where the number of variables
is part of the indexing (compositions have length `n`, permutations have their last descent at most
`n`, and the elementary basis is available).

EXAMPLES::

    sage: from schubmult.sage import PolynomialAlgebra
    sage: A = PolynomialAlgebra(QQ); A
    Polynomial ring in x0, x1, ... over Rational Field with combinatorial bases
    sage: S = A.schubert(); k = A.key(); x = A.monomial()
    sage: S[3, 1, 2]
    S[3, 1, 2]
    sage: k(S[3, 1, 2])
    k[2]
    sage: x(k[2, 0, 1])
    x[2, 0, 1] + x[2, 1]
    sage: k[2, 0, 1].expand()
    x0^2*x1 + x0^2*x2
    sage: S[2, 1, 4, 3] * S[1, 3, 2]
    S[2, 3, 4, 1] + S[2, 4, 1, 3] + S[3, 1, 4, 2] + S[4, 1, 2, 3]
    sage: A.grothendieck()(S[1, 3, 2])
    G[1, 3, 2] - G[2, 3, 1]

The bases coerce into one another, and Sage's own Schubert and key polynomials coerce in::

    sage: k[1, 0, 2] + S[1, 3, 2]
    k[0, 1] + k[1, 0, 2]
    sage: A.forest()(SchubertPolynomialRing(QQ)([3, 1, 2]))
    P[2]
    sage: A.fundamental_slide()(KeyPolynomials(QQ)([1, 0, 2]))
    F[1, 0, 2] + F[2, 0, 1]

In a fixed number of variables::

    sage: B = PolynomialAlgebra(QQ, 3); B
    Polynomial ring in x0, x1, x2 over Rational Field with combinatorial bases
    sage: E = B.elementary(); E[1, 0, 2]
    E[1, 0, 2]
    sage: E[1, 0, 2].expand()
    x0^2*x1 + x0^2*x2 + x0*x1*x2
    sage: B.key()(E[1, 0, 2])
    k[1, 1, 1] + k[2, 0, 1]
    sage: B.key()[2, 0, 1] * B.key()[0, 1, 0]
    k[2, 1, 1] + k[2, 2, 0] + k[3, 0, 1]
"""

from sage.categories.algebras import Algebras
from sage.categories.realizations import Category_realization_of_parent
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.integer_vector import IntegerVectors
from sage.combinat.key_polynomial import OperatorPolynomial, OperatorPolynomialBasis
from sage.combinat.permutation import Permutation, Permutations
from sage.combinat.schubert_polynomial import SchubertPolynomialRing_xbasis
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.infinite_polynomial_element import InfinitePolynomial
from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

from ._common import X_LETTER, genset, to_sage_perm, to_schubmult_perm
from ._convert import parse_sage_name, symengine_to_sage


def _trim(alpha):
    alpha = list(alpha)
    while alpha and alpha[-1] == 0:
        alpha.pop()
    return tuple(alpha)


class _Backend:
    """A schubmult ``PolynomialBasis`` with its conversions, on ``{schubmult key: SymEngine coefficient}``
    dicts whose keys all refer to the same number of variables."""

    def __init__(self, make_basis):
        from schubmult.rings.polynomial_algebra import MonomialBasis, PolynomialAlgebra

        self.basis = make_basis(genset(X_LETTER))
        self.ring = PolynomialAlgebra(self.basis)
        self.monomial = MonomialBasis(genset(X_LETTER))

    def expand(self, key):
        from schubmult.symbolic import S

        return self.ring.from_dict({key: S.One}).expand()

    def to_monomials(self, dct):
        return self.basis.transition(self.monomial)(dct)

    def from_monomials(self, mono):
        return self.monomial.transition(self.basis)(mono)

    def product(self, a, b):
        return self.basis.product(a, b)


def _bases():
    from schubmult.rings.polynomial_algebra import (
        ElemSymPolyBasis,
        ForestPolyBasis,
        FundamentalSlidePolyBasis,
        GlidePolyBasis,
        GrothendieckPolyBasis,
        GrovePolyBasis,
        KeyPolyBasis,
        LascouxPolyBasis,
        MonomialBasis,
        MonomialSlidePolyBasis,
        SchubertPolyBasis,
    )

    return {
        # name: (realization class, prefix, description, schubmult basis factory)
        "monomial": (_CompositionBasis, "x", "monomial", MonomialBasis),
        "schubert": (_PermutationBasis, "S", "Schubert", SchubertPolyBasis),
        "grothendieck": (_PermutationBasis, "G", "Grothendieck (beta = 1)", GrothendieckPolyBasis),
        "key": (_CompositionBasis, "k", "key polynomial", KeyPolyBasis),
        "fundamental_slide": (_CompositionBasis, "F", "fundamental slide polynomial", FundamentalSlidePolyBasis),
        "monomial_slide": (_CompositionBasis, "M", "monomial slide polynomial", MonomialSlidePolyBasis),
        "forest": (_CompositionBasis, "P", "forest polynomial", ForestPolyBasis),
        "glide": (_CompositionBasis, "Gl", "glide polynomial (beta = 1)", GlidePolyBasis),
        "lascoux": (_CompositionBasis, "L", "Lascoux polynomial (beta = 1)", LascouxPolyBasis),
        "grove": (_CompositionBasis, "Gr", "grove polynomial (beta = 1)", lambda gs: GrovePolyBasis(gs, beta=1)),
        "elementary": (_ElementaryBasis, "E", "elementary symmetric", ElemSymPolyBasis),
    }


class PolynomialAlgebra(UniqueRepresentation, Parent):
    r"""
    The polynomial ring `R[x_0, x_1, \ldots]`, or `R[x_0, \ldots, x_{n-1}]` if ``n`` is given, with its
    combinatorial bases as realizations.

    EXAMPLES::

        sage: from schubmult.sage import PolynomialAlgebra
        sage: A = PolynomialAlgebra(ZZ); A
        Polynomial ring in x0, x1, ... over Integer Ring with combinatorial bases
        sage: sorted(A.basis_names())
        ['forest', 'fundamental_slide', 'glide', 'grothendieck', 'grove', 'key', 'lascoux', 'monomial', 'monomial_slide', 'schubert']
        sage: A.basis('lascoux')
        Polynomial ring in x0, x1, ... over Integer Ring in the Lascoux polynomial (beta = 1) basis
        sage: A.elementary()
        Traceback (most recent call last):
        ...
        ValueError: the elementary symmetric basis depends on the number of variables: use PolynomialAlgebra(R, n)
        sage: TestSuite(A).run()

        sage: B = PolynomialAlgebra(ZZ, 2); B
        Polynomial ring in x0, x1 over Integer Ring with combinatorial bases
        sage: B.number_of_variables()
        2
        sage: 'elementary' in B.basis_names()
        True
        sage: TestSuite(B).run()
    """

    @staticmethod
    def __classcall_private__(cls, R, n=None):  # noqa: PLW0211  (Sage's classcall protocol)
        if n is not None:
            if n not in ZZ or n < 1:
                raise ValueError("the number of variables must be a positive integer")
            n = Integer(n)
        return super().__classcall__(cls, R, n)

    def __init__(self, R, n=None):
        self._n = None if n is None else int(n)
        Parent.__init__(self, base=R, category=Algebras(R).Commutative().WithRealizations())

    def _repr_(self):
        return f"{self._ring_name()} with combinatorial bases"

    def _ring_name(self):
        variables = "x0, x1, ..." if self._n is None else ", ".join(f"{X_LETTER}{i}" for i in range(self._n))
        return f"Polynomial ring in {variables} over {self.base_ring()}"

    def number_of_variables(self):
        """The number of variables, or ``None`` for the ring in infinitely many variables."""
        return self._n

    def a_realization(self):
        return self.monomial()

    def basis_names(self):
        """The names of the bases available in this ring."""
        return [name for name in _bases() if name != "elementary" or self._n is not None]

    @cached_method
    def basis(self, name):
        """The realization called ``name`` (see :meth:`basis_names`)."""
        if name == "elementary" and self._n is None:
            raise ValueError("the elementary symmetric basis depends on the number of variables: use PolynomialAlgebra(R, n)")
        if name not in _bases():
            raise ValueError(f"unknown basis {name!r}; choose from {sorted(_bases())}")
        return _bases()[name][0](self, name)

    def monomial(self):
        r"""
        The monomial basis `x^\alpha`.

        EXAMPLES::

            sage: from schubmult.sage import PolynomialAlgebra
            sage: x = PolynomialAlgebra(QQ).monomial()
            sage: x[2, 0, 1] * x[0, 1]
            x[2, 1, 1]
            sage: x[2, 0, 1].expand()
            x0^2*x2
        """
        return self.basis("monomial")

    def schubert(self):
        r"""
        Schubert polynomials `\mathfrak S_w`, indexed by permutations.

        EXAMPLES::

            sage: from schubmult.sage import PolynomialAlgebra
            sage: S = PolynomialAlgebra(QQ).schubert()
            sage: S(Permutation([3, 1, 2])) == S[3, 1, 2] == S([3, 1, 2, 4])
            True
            sage: S[3, 1, 2] * S[1, 3, 2]
            S[3, 2, 1] + S[4, 1, 2, 3]
            sage: S[2, 3, 1].expand()
            x0*x1
            sage: S[1, 3, 2].degree()
            1

        In `n` variables the last descent must be at most `n`::

            sage: PolynomialAlgebra(QQ, 2).schubert()([1, 2, 4, 3])
            Traceback (most recent call last):
            ...
            ValueError: [1, 2, 4, 3] has a descent beyond the 2 variables of Polynomial ring in x0, x1 over Rational Field with combinatorial bases
        """
        return self.basis("schubert")

    def grothendieck(self):
        r"""
        Grothendieck polynomials `\mathfrak G_w` at `\beta = 1`, indexed by permutations.

        EXAMPLES::

            sage: from schubmult.sage import PolynomialAlgebra
            sage: A = PolynomialAlgebra(QQ); G = A.grothendieck()
            sage: G[1, 3, 2].expand()
            x0*x1 + x0 + x1
            sage: A.schubert()(G[1, 3, 2])
            S[1, 3, 2] + S[2, 3, 1]
            sage: G[2, 1] * G[2, 1]
            G[3, 1, 2]
        """
        return self.basis("grothendieck")

    def key(self):
        r"""
        Key polynomials (Demazure characters) `\kappa_\alpha`.

        EXAMPLES::

            sage: from schubmult.sage import PolynomialAlgebra
            sage: k = PolynomialAlgebra(QQ).key()
            sage: k[0, 2].expand()
            x0^2 + x0*x1 + x1^2
            sage: k(KeyPolynomials(QQ)([0, 2]))
            k[0, 2]
        """
        return self.basis("key")

    def fundamental_slide(self):
        r"""
        Fundamental slide polynomials `\mathfrak F_\alpha` (Assaf-Searles).

        EXAMPLES::

            sage: from schubmult.sage import PolynomialAlgebra
            sage: A = PolynomialAlgebra(QQ); F = A.fundamental_slide()
            sage: F(A.schubert()[2, 1, 5, 3, 4])
            F[1, 0, 2] + F[2, 0, 1] + F[3]
        """
        return self.basis("fundamental_slide")

    def monomial_slide(self):
        r"""
        Monomial slide polynomials `\mathfrak M_\alpha` (Assaf-Searles).

        EXAMPLES::

            sage: from schubmult.sage import PolynomialAlgebra
            sage: A = PolynomialAlgebra(QQ); M = A.monomial_slide()
            sage: M(A.fundamental_slide()[1, 0, 2])
            M[1, 0, 2] + M[1, 1, 1]
        """
        return self.basis("monomial_slide")

    def forest(self):
        r"""
        Forest polynomials `\mathfrak P_F` (Nadeau-Spink-Tewari), indexed by weak compositions.

        EXAMPLES::

            sage: from schubmult.sage import PolynomialAlgebra
            sage: A = PolynomialAlgebra(QQ); P = A.forest()
            sage: P(A.schubert()[2, 1, 5, 3, 4])
            P[1, 0, 2]
        """
        return self.basis("forest")

    def glide(self):
        r"""
        Glide polynomials (Pechenik-Searles) at `\beta = 1`, the K-theoretic fundamental slides.

        EXAMPLES::

            sage: from schubmult.sage import PolynomialAlgebra
            sage: A = PolynomialAlgebra(QQ); Gl = A.glide()
            sage: Gl[0, 2].expand()
            x0^2*x1 + x0*x1^2 + x0^2 + x0*x1 + x1^2
        """
        return self.basis("glide")

    def lascoux(self):
        r"""
        Lascoux polynomials at `\beta = 1`, the K-theoretic key polynomials.

        EXAMPLES::

            sage: from schubmult.sage import PolynomialAlgebra
            sage: A = PolynomialAlgebra(QQ); L = A.lascoux()
            sage: A.glide()(L[0, 2])
            Gl[0, 2]
        """
        return self.basis("lascoux")

    def grove(self):
        r"""
        Grove polynomials at `\beta = 1`, the K-theoretic forest polynomials.

        EXAMPLES::

            sage: from schubmult.sage import PolynomialAlgebra
            sage: A = PolynomialAlgebra(QQ); Gr = A.grove()
            sage: Gr[0, 2].expand()
            x0^2*x1 + x0*x1^2 + x0^2 + x0*x1 + x1^2
        """
        return self.basis("grove")

    def elementary(self):
        r"""
        Products of elementary symmetric polynomials, a basis of the ring in `n` variables.

        The index is a tuple `(a_1, \ldots, a_{n-1}, b_1, \ldots, b_r)`, `r \geq 1`, standing for
        `\prod_{j=1}^{n-1} e_{a_j}(x_0, \ldots, x_{j-1}) \cdot \prod_i e_{b_i}(x_0, \ldots, x_{n-1})` with
        `b_1 \leq \cdots \leq b_r` (a single `b_1 = 0` when there is no full-alphabet factor): the
        basis `\{\prod_{j<n} e_{a_j}(x_0..x_{j-1})\}` of the polynomial ring as a free module over the
        symmetric polynomials, times the monomials in `e_1, \ldots, e_n`. Since the meaning of an index
        depends on `n`, this basis exists only in ``PolynomialAlgebra(R, n)``.

        EXAMPLES::

            sage: from schubmult.sage import PolynomialAlgebra
            sage: B = PolynomialAlgebra(QQ, 3); E = B.elementary()
            sage: E[1, 0, 2].expand()                 # e_1(x0) e_0(x0, x1) e_2(x0, x1, x2)
            x0^2*x1 + x0^2*x2 + x0*x1*x2
            sage: E[0, 0, 1, 1].expand()              # e_1(x0, x1, x2)^2
            x0^2 + 2*x0*x1 + x1^2 + 2*x0*x2 + 2*x1*x2 + x2^2
            sage: E([1, 0, 2, 0, 1]) == E[1, 0, 1, 2]  # e_0 factors drop out, full-alphabet degrees are sorted
            True
            sage: E(B.monomial()[2])
            -E[0, 2, 0] + E[1, 1, 0]
            sage: B.schubert()(E[0, 0, 1])
            S[1, 2, 4, 3]
        """
        return self.basis("elementary")

    class Bases(Category_realization_of_parent):
        def super_categories(self):
            A = self.base()
            return [A.Realizations(), Algebras(A.base_ring()).Commutative().WithBasis().Filtered().Realizations()]


class PolynomialAlgebraBasis(CombinatorialFreeModule):
    """
    A basis of :func:`PolynomialAlgebra`; subclasses fix the index set.

    EXAMPLES::

        sage: from schubmult.sage import PolynomialAlgebra
        sage: k = PolynomialAlgebra(QQ).key(); k
        Polynomial ring in x0, x1, ... over Rational Field in the key polynomial basis
        sage: k([2, 0, 1, 0, 0])
        k[2, 0, 1]
        sage: k.one()
        k[]
        sage: TestSuite(k).run()
        sage: TestSuite(PolynomialAlgebra(QQ, 2).schubert()).run()
    """

    def __init__(self, A, name):
        self._name = name
        _, self._prefix, self._description, make = _bases()[name]
        self._backend = _Backend(make)
        CombinatorialFreeModule.__init__(self, A.base_ring(), self._index_set(A), prefix=self._prefix, bracket=False, category=A.Bases(), sorting_key=self._sorting_key)

    def _repr_(self):
        return f"{self.realization_of()._ring_name()} in the {self._description} basis"

    def _realization_name(self):
        return self._name

    def _n(self):
        return self.realization_of()._n

    # ---- indexing (subclasses) -------------------------------------------------------------

    def _index_set(self, A):
        raise NotImplementedError

    @staticmethod
    def _sorting_key(key):
        raise NotImplementedError

    def _key(self, x):
        """User input (list, Sage index element, ...) -> index of this basis."""
        raise NotImplementedError

    def _variables_of(self, key):
        """Number of variables the basis element indexed by ``key`` needs."""
        raise NotImplementedError

    def _schubmult_key(self, key, N):
        """The schubmult key for ``key`` in ``N`` variables."""
        raise NotImplementedError

    def _from_schubmult_key(self, key):
        raise NotImplementedError

    def _some_keys(self):
        raise NotImplementedError

    # ---- conversions -----------------------------------------------------------------------

    def _variables(self, keys):
        """The number of variables to run a schubmult computation on ``keys`` in."""
        if self._n() is not None:
            return self._n()
        return max((self._variables_of(k) for k in keys), default=0)

    def _scalar(self, c):
        import symengine

        c = symengine.sympify(c)
        R = self.base_ring()
        return R(int(c)) if c.is_Integer else R(QQ((int(c.p), int(c.q))))

    def _from_schubmult(self, dct):
        out = {}
        for k, c in dct.items():
            if c == 0:
                continue
            key = self._from_schubmult_key(k)
            out[key] = out[key] + self._scalar(c) if key in out else self._scalar(c)
        return self._from_dict(out, remove_zeros=True)

    def _to_schubmult(self, elem, N=None):
        from schubmult.symbolic import sympify

        N = self._variables(elem.support()) if N is None else N
        return {self._schubmult_key(k, N): sympify(str(c)) for k, c in elem}

    def _to_monomials(self, elem):
        """``{exponents (all of one length): coefficient}`` of an element of this basis."""
        return self._backend.to_monomials(self._to_schubmult(elem))

    def _from_monomials(self, mono):
        mono = self._pad_monomials(mono)
        return self._from_schubmult(self._backend.from_monomials(mono))

    def _pad_monomials(self, mono):
        n = self._n()
        if n is None:
            n = max((len(_trim(k)) for k in mono), default=0)
        elif any(len(_trim(k)) > n for k in mono):
            raise ValueError(f"a polynomial in more than the {n} variables of {self.realization_of()}")
        return {_trim(k) + (0,) * (n - len(_trim(k))): c for k, c in mono.items()}

    def _from_polynomial(self, p):
        """A Sage polynomial, its indexed variables read as ``x0, x1, ...`` by index."""
        from schubmult.symbolic import sympify

        if isinstance(p, InfinitePolynomial):
            p = p.polynomial()
        positions = [parse_sage_name(v)[1] for v in p.parent().variable_names()]
        mono = {}
        for exps, c in p.dict().items():
            key = [0] * (1 + max((i for i, e in zip(positions, exps) if e), default=-1))
            for i, e in zip(positions, exps):
                if e:
                    key[i] = e
            key = tuple(key)
            mono[key] = mono.get(key, 0) + sympify(str(c))
        return self._from_monomials(mono)

    # ---- algebra structure -----------------------------------------------------------------

    @cached_method
    def one_basis(self):
        return self._key([])

    def product_on_basis(self, left, right):
        N = self._variables([left, right])
        return self._from_schubmult(self._backend.product(self._schubmult_key(left, N), self._schubmult_key(right, N)))

    def __getitem__(self, key):
        """``k[2, 0, 1]``, ``S[3, 1, 2]``: the basis element with that index."""
        if isinstance(key, int | Integer):
            key = [key]
        return self.monomial(self._key(list(key)))

    def _element_constructor_(self, x):
        if isinstance(x, list | tuple | Permutation) or (hasattr(x, "parent") and x.parent() is self._indices):
            return self.monomial(self._key(x))
        if isinstance(x, Polynomial):
            x = PolynomialRing(x.base_ring(), 1, x.parent().variable_name())(x)
        if isinstance(x, OperatorPolynomial):  # Sage key/atom polynomials: their variables are the x's
            x = x.expand()
        if isinstance(x, MPolynomial | InfinitePolynomial):
            return self._from_polynomial(x)
        parent = getattr(x, "parent", None)
        if parent is not None:
            parent = parent()
            if isinstance(parent, PolynomialAlgebraBasis) and parent.realization_of() is self.realization_of():
                return self._from_monomials(parent._to_monomials(x))
            if isinstance(parent, PolynomialAlgebraBasis | SchubertPolynomialRing_xbasis):
                return self._from_polynomial(x.expand())
            if parent is self.base_ring() or self.base_ring().has_coerce_map_from(parent):
                return self.from_base_ring(self.base_ring()(x))
        raise TypeError(f"do not know how to make an element of {self} from {x!r}")

    def _coerce_map_from_(self, S):
        if isinstance(S, PolynomialAlgebraBasis):
            same = S.realization_of() is self.realization_of()
            fits = S._n() is not None and (self._n() is None or S._n() <= self._n())  # fewer variables embed
            return same or (fits and self.base_ring().has_coerce_map_from(S.base_ring()))
        if isinstance(S, SchubertPolynomialRing_xbasis | OperatorPolynomialBasis):
            return self._n() is None and self.base_ring().has_coerce_map_from(S.base_ring())
        return super()._coerce_map_from_(S)

    def some_elements(self):
        keys = self._some_keys()
        return [self.one(), self.one() + 2 * self(keys[0]), self(keys[1]) - self(keys[2])]

    class Element(CombinatorialFreeModule.Element):
        def expand(self):
            r"""
            The polynomial in ``x0, x1, ...``: in the `n` variables of the ring, or in as many as the
            support needs (at least one).

            EXAMPLES::

                sage: from schubmult.sage import PolynomialAlgebra
                sage: A = PolynomialAlgebra(QQ)
                sage: A.key()[1, 0, 2].expand()
                x0^2*x1 + x0*x1^2 + x0^2*x2 + x0*x1*x2 + x0*x2^2
                sage: A.lascoux()[0, 2].expand()
                x0^2*x1 + x0*x1^2 + x0^2 + x0*x1 + x1^2
                sage: A.monomial().one().expand().parent()
                Multivariate Polynomial Ring in x0 over Rational Field
                sage: PolynomialAlgebra(QQ, 3).monomial().one().expand().parent()
                Multivariate Polynomial Ring in x0, x1, x2 over Rational Field
            """
            P = self.parent()
            N = P._variables(self.support())
            T = PolynomialRing(P.base_ring(), max(N, 1), [f"{X_LETTER}{i}" for i in range(max(N, 1))])
            gens = T.gens()

            def variable(letter, i):  # noqa: ARG001
                return gens[i - 1]

            def scalar(q):
                return T(q) if isinstance(q, int) else T(QQ(q.numerator) / QQ(q.denominator))

            return sum((c * symengine_to_sage(P._backend.expand(P._schubmult_key(k, N)), variable, scalar) for k, c in self), T.zero())


class _CompositionBasis(PolynomialAlgebraBasis):
    """Indexed by weak compositions: trailing zeros dropped, or padded to the `n` variables of the ring."""

    def _index_set(self, A):
        return IntegerVectors() if A._n is None else IntegerVectors(k=A._n)

    @staticmethod
    def _sorting_key(key):
        return (sum(key), list(key))

    def _key(self, x):
        if isinstance(x, Permutation):
            raise TypeError(f"permutations index the Schubert and Grothendieck bases, not the {self._description} basis")
        alpha = _trim(x)
        n = self._n()
        if n is None:
            return self._indices(list(alpha))
        if len(alpha) > n:
            raise ValueError(f"{list(x)} has more than the {n} variables of {self.realization_of()}")
        return self._indices([*alpha, *[0] * (n - len(alpha))])

    def _variables_of(self, key):
        return len(_trim(key))

    def _schubmult_key(self, key, N):
        alpha = _trim(key)
        return alpha + (0,) * (N - len(alpha))

    def _from_schubmult_key(self, key):
        return self._key(key)

    def degree_on_basis(self, key):
        return sum(key)

    def _some_keys(self):
        n = self._n()
        return [[1], [2, 0, 1], [0, 1]] if n is None or n >= 3 else ([[1], [2, 0], [0, 1]] if n == 2 else [[1], [2], [3]])


class _PermutationBasis(PolynomialAlgebraBasis):
    """Indexed by permutations; the schubmult key is ``(w, N)`` for the number of variables ``N`` in play."""

    def _index_set(self, A):  # noqa: ARG002
        return Permutations()

    @staticmethod
    def _sorting_key(key):
        return (len(key), list(key))

    def _key(self, x):
        w = to_sage_perm(x)
        n = self._n()
        if n is not None and self._variables_of(w) > n:
            raise ValueError(f"{list(x)} has a descent beyond the {n} variables of {self.realization_of()}")
        return w

    def _variables_of(self, key):
        return max(key.descents(), default=0)  # position of the last descent, 1-based

    def _schubmult_key(self, key, N):
        return (to_schubmult_perm(key), N)

    def _from_schubmult_key(self, key):
        return to_sage_perm(key[0])

    def degree_on_basis(self, key):
        return key.length()

    def _some_keys(self):
        n = self._n()
        return [[2, 1], [3, 1, 2], [1, 3, 2]] if n is None or n >= 2 else [[2, 1], [3, 1, 2], [4, 1, 2, 3]]


class _ElementaryBasis(PolynomialAlgebraBasis):
    """Products of elementary symmetric polynomials in exactly the ``n`` variables of the ring; see
    :meth:`PolynomialAlgebra.elementary` for the meaning of an index."""

    def _index_set(self, A):  # noqa: ARG002
        return IntegerVectors()

    @staticmethod
    def _sorting_key(key):
        return (sum(key), list(key))

    def _canonical(self, tup):
        n = self._n()
        tup = list(tup)
        lower = tup[: n - 1] + [0] * (n - 1 - len(tup[: n - 1]))
        if any(a > j for j, a in enumerate(lower, start=1)):
            raise ValueError(f"{tup}: e_a(x0..x_(j-1)) needs a <= j")
        top = sorted(a for a in tup[n - 1 :] if a)
        if top and top[-1] > n:
            raise ValueError(f"{tup}: e_a in {n} variables needs a <= {n}")
        return tuple(lower + (top or [0]))

    def _key(self, x):
        return self._indices(list(self._canonical(x)))

    def _variables_of(self, key):  # noqa: ARG002
        return self._n()

    def _schubmult_key(self, key, N):  # noqa: ARG002
        return (tuple(key), self._n())

    def _from_schubmult_key(self, key):
        return self._key(key[0])

    def degree_on_basis(self, key):
        return sum(key)

    def product_on_basis(self, left, right):
        # schubmult multiplies this basis through the monomial one; do the same, in the ring's n variables
        a = self._to_monomials(self.monomial(left))
        b = self._to_monomials(self.monomial(right))
        mono = {}
        for ka, ca in a.items():
            for kb, cb in b.items():
                k = tuple(i + j for i, j in zip(ka, kb))
                mono[k] = mono.get(k, 0) + ca * cb
        return self._from_monomials(mono)

    def _some_keys(self):
        n = self._n()
        return [[0] * (n - 1) + [1], [1] + [0] * (n - 1), [0] * (n - 1) + [1, 1]]
