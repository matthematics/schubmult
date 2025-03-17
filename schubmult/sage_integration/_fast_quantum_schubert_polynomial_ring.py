import schubmult.schubmult_q as sq


# sage.doctest: needs sage.combinat sage.modules
r"""
QuantumSchubert Polynomials


See :wikipedia:`QuantumSchubert_polynomial` and
`SymmetricFunctions.com <https://www.symmetricfunctions.com/schubert.htm#schubert>`_.
QuantumSchubert polynomials are representatives of cohomology classes in flag varieties.
In `n` variables, they are indexed by permutations `w \in S_n`. They also form
a basis for the coinvariant ring of the `S_n` action on
`\ZZ[x_1, x_2, \ldots, x_n]`.

EXAMPLES::

	sage: X = FastQuantumSchubertPolynomialRing(ZZ)
	sage: w = [1,2,5,4,3];  # a list representing an element of `S_5`
	sage: X(w)
	X[1, 2, 5, 4, 3]

This can be expanded in terms of polynomial variables::

	sage: X(w).expand()
	x0^2*x1 + x0*x1^2 + x0^2*x2 + 2*x0*x1*x2 + x1^2*x2
	 + x0*x2^2 + x1*x2^2 + x0^2*x3 + x0*x1*x3 + x1^2*x3
	 + x0*x2*x3 + x1*x2*x3 + x2^2*x3

We can also convert back from polynomial variables. For example,
the longest permutation is a single term. In `S_5`, this is the
element (in one line notation) `w_0 = 54321`::

	sage: w0 = [5,4,3,2,1]
	sage: R.<x0, x1, x2, x3, x4> = PolynomialRing(ZZ)
	sage: Sw0 = X(x0^4*x1^3*x2^2*x3);  Sw0
	X[5, 4, 3, 2, 1]

The polynomials also have the property that if the indexing permutation `w` is
multiplied by a simple transposition `s_i = (i, i+1)` such that the length of
`w` is more than the length of `ws_i`, then the FastQuantumSchubert
polynomial of the permutation `ws_i` is computed by applying the divided
difference operator :meth:`~QuantumSchubertPolynomial_class.divided_difference` to
the polynomial indexed by `w`. For example, applying the divided difference
operator `\partial_2` to the FastQuantumSchubert polynomial `\mathfrak{S}_{w_0}`::

	sage: Sw0.divided_difference(2)
	X[5, 3, 4, 2, 1]

We can also check the properties listed in :wikipedia:`QuantumSchubert_polynomial`::

	sage: X([1,2,3,4,5])  # the identity in one-line notation
	X[1]
	sage: X([1,3,2,4,5]).expand()  # the transposition swapping 2 and 3
	x0 + x1
	sage: X([2,4,5,3,1]).expand()
	x0^2*x1^2*x2*x3 + x0^2*x1*x2^2*x3 + x0*x1^2*x2^2*x3

	sage: w = [4,5,1,2,3]
	sage: s = SymmetricFunctions(QQ).schur()
	sage: s[3,3].expand(2)
	x0^3*x1^3
	sage: X(w).expand()
	x0^3*x1^3
"""
# ****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.combinat.free_module import CombinatorialFreeModule

# from sage.combinat.key_polynomial import KeyPolynomial
from sage.combinat.permutation import Permutations, Permutation, from_lehmer_code
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.infinite_polynomial_element import InfinitePolynomial
from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

lazy_import("sage.libs.symmetrica", "all", as_="symmetrica")


# def FastQuantumSchubertPolynomialRing(R):
# """
# Return the FastQuantumSchubert polynomial ring over ``R`` on the X basis.
#
# This is the basis made of the FastQuantumSchubert polynomials.
#
# EXAMPLES::
#
# sage: X = FastQuantumSchubertPolynomialRing(ZZ); X
# QuantumSchubert polynomial ring with X basis over Integer Ring
# sage: TestSuite(X).run()
# sage: X(1)
# X[1]
# sage: X([1,2,3])*X([2,1,3])
# X[2, 1]
# sage: X([2,1,3])*X([2,1,3])
# X[3, 1, 2]
# sage: X([2,1,3])+X([3,1,2,4])
# X[2, 1] + X[3, 1, 2]
# sage: a = X([2,1,3])+X([3,1,2,4])
# sage: a^2
# X[3, 1, 2] + 2*X[4, 1, 2, 3] + X[5, 1, 2, 3, 4]
# """
# return FastQuantumSchubertPolynomialRing_xbasis(R)
#
def FastQuantumSchubertPolynomialRing(R):
    QR = PolynomialRing(R, 100, "q_")
    return FastQuantumSchubertPolynomialRing_xbasis(QR)


class FastQuantumSchubertPolynomial_class(CombinatorialFreeModule.Element):
    def expand(self):
        """
        EXAMPLES::

                sage: X = FastQuantumSchubertPolynomialRing(ZZ)
                sage: X([2,1,3]).expand()
                x0
                sage: [X(p).expand() for p in Permutations(3)]
                [1, x0 + x1, x0, x0*x1, x0^2, x0^2*x1]

        TESTS:

        Calling .expand() should always return an element of an
        MPolynomialRing::

                sage: X = FastQuantumSchubertPolynomialRing(ZZ)
                sage: f = X([1]); f
                X[1]
                sage: type(f.expand())
                <class Fast'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular'>
                sage: f.expand()
                1
                sage: f = X([1,2])
                sage: type(f.expand())
                <class Fast'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular'>
                sage: f = X([1,3,2,4])
                sage: type(f.expand())
                <class Fast'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular'>

        Now we check for correct handling of the empty
        permutation (:issue:`23443`)::

                sage: X([1]).expand() * X([2,1]).expand()
                x0
        """
        p = symmetrica.t_SCHUBERT_POLYNOM(self)
        if not isinstance(p, MPolynomial):
            R = PolynomialRing(self.parent().base_ring(), 1, "x0")
            p = R(p)
        return p


class FastQuantumSchubertPolynomialRing_xbasis(CombinatorialFreeModule):
    Element = FastQuantumSchubertPolynomial_class

    def __init__(self, R):
        """
        EXAMPLES::

                sage: X = FastQuantumSchubertPolynomialRing(QQ)
                sage: X == loads(dumps(X))
                True
        """
        self._name = "QuantumSchubert polynomial ring with X basis"
        self._repr_option_bracket = False
        cat = GradedAlgebrasWithBasis(R)  # CoalgebrasWithBasis(R).Graded()
        CombinatorialFreeModule.__init__(
            self, R, Permutations(), category=cat, prefix="X"
        )

    @cached_method
    def one_basis(self):
        """
        Return the index of the unit of this algebra.

        EXAMPLES::

                sage: X = FastQuantumSchubertPolynomialRing(QQ)
                sage: X.one()  # indirect doctest
                X[1]
        """
        return self._indices([1])

    def _element_constructor_(self, x):
        """
        Coerce x into ``self``.

        EXAMPLES::

                sage: X = FastQuantumSchubertPolynomialRing(QQ)
                sage: X._element_constructor_([2,1,3])
                X[2, 1]
                sage: X._element_constructor_(Permutation([2,1,3]))
                X[2, 1]

                sage: R.<x1, x2, x3> = QQ[]
                sage: X(x1^2*x2)
                X[3, 2, 1]

                sage: S.<x> = InfinitePolynomialRing(QQ)
                sage: X(x[0]^2*x[1])
                X[3, 2, 1]
                sage: X(x[0]*x[1]^2*x[2]^2*x[3] + x[0]^2*x[1]^2*x[2]*x[3] + x[0]^2*x[1]*x[2]^2*x[3])
                X[2, 4, 5, 3, 1]

                sage: from sage.combinat.key_polynomial import KeyPolynomialBasis
                sage: k = KeyPolynomialBasis(QQ)
                sage: X(k([3,2,1]))
                X[4, 3, 2, 1]

        TESTS:

        We check that :issue:`12924` is fixed::

                sage: X = FastQuantumSchubertPolynomialRing(QQ)
                sage: X._element_constructor_([1,2,1])
                Traceback (most recent call last):
                ...
                ValueError: the input [1, 2, 1] is not a valid permutation

        Now we check for correct handling of the empty
        permutation (:issue:`23443`)::

                sage: X([])
                X[1]

        Check the round trip from key polynomials::

                sage: k = KeyPolynomials(ZZ)
                sage: X = FastQuantumSchubertPolynomialRing(ZZ)
                sage: it = iter(Permutations())
                sage: for _ in range(50):
                ....:     P = next(it)
                ....:     assert X(k(X(P))) == X(P), P
        """
        if isinstance(x, list):
            # checking the input to avoid symmetrica crashing Sage, see trac 12924
            if x not in Permutations():
                raise ValueError(f"the input {x} is not a valid permutation")
            perm = Permutation(x).remove_extra_fixed_points()
            return self._from_dict({perm: self.base_ring().one()})
        elif isinstance(x, Permutation):
            perm = x.remove_extra_fixed_points()
            return self._from_dict({perm: self.base_ring().one()})
        elif isinstance(x, MPolynomial):
            return symmetrica.t_POLYNOM_SCHUBERT(x)
        elif isinstance(x, InfinitePolynomial):
            R = x.polynomial().parent()
            # massage the term order to be what symmetrica expects
            S = PolynomialRing(R.base_ring(), names=list(map(repr, reversed(R.gens()))))
            return symmetrica.t_POLYNOM_SCHUBERT(S(x.polynomial()))
        # elif isinstance(x, KeyPolynomial):
        # return self(x.expand())
        else:
            raise TypeError

    def some_elements(self):
        """
        Return some elements.

        EXAMPLES::

                sage: X = FastQuantumSchubertPolynomialRing(QQ)
                sage: X.some_elements()
                [X[1], X[1] + 2*X[2, 1], -X[3, 2, 1] + X[4, 2, 1, 3]]
        """
        return [
            self.one(),
            self([1]) + 2 * self([2, 1]),
            self([4, 2, 1, 3]) - self([3, 2, 1]),
        ]

    def product_on_basis(self, left, right):
        """
        EXAMPLES::

                sage: p1 = Permutation([3,2,1])
                sage: p2 = Permutation([2,1,3])
                sage: X = FastQuantumSchubertPolynomialRing(QQ)
                sage: X.product_on_basis(p1,p2)
                X[4, 2, 1, 3]
        """

        # return symmetrica.mult_schubert_schubert(left, right)
        return sum(
            [
                self.base_ring()(v) * self(Permutation(list(k)))
                for k, v in sq.schubmult_db(
                    {tuple(left): 1}, tuple(right), list(self.base_ring().gens())
                ).items()
            ]
        )
