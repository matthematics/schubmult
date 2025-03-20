import schubmult.schubmult_q as sq
import schubmult.schubmult_q_yz as qyz
from sympy import sympify
import symengine as syme


from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.combinat.free_module import CombinatorialFreeModule

from sage.combinat.permutation import Permutations, Permutation
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import
from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_base
from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

lazy_import("sage.libs.symmetrica", "all", as_="symmetrica")


def FastQuantumSchubertPolynomialRing(R, num_vars, varname, q_varname="q_"):
    QR = PolynomialRing(R, num_vars, q_varname)
    return FastQuantumSchubertPolynomialRing_xbasis(QR, num_vars, varname, q_varname)


class FastQuantumSchubertPolynomial_class(CombinatorialFreeModule.Element):
    def expand(self):
        return sum(
            [
                self.parent()._polynomial_ring(
                    qyz.schubpoly_quantum(
                        tuple(k),
                        self.parent()._polynomial_ring.gens(),
                        [0 for i in range(100)],
                        self.parent()._q_ring.gens(),
                        v,
                    )
                )
                for k, v in self.monomial_coefficients().items()
            ]
        )


class FastQuantumSchubertPolynomialRing_xbasis(CombinatorialFreeModule):
    Element = FastQuantumSchubertPolynomial_class

    def __init__(self, R, num_vars, varname, q_varname):
        """
        EXAMPLES::

                sage: X = FastQuantumSchubertPolynomialRing(QQ)
                sage: X == loads(dumps(X))
                True
        """
        self._name = "QuantumSchubert polynomial ring with X basis"
        self._repr_option_bracket = False
        cat = GradedAlgebrasWithBasis(R).Commutative()  # CoalgebrasWithBasis(R).Graded()
        CombinatorialFreeModule.__init__(
            self, R, Permutations(), category=cat, prefix=f"S^{q_varname}({varname})"
        )
        self._q_ring = R
        self._base_varname = varname
        self._q_varname = q_varname
        self._polynomial_ring = PolynomialRing(R, num_vars, varname)
        self._populate_coercion_lists_()

    def _coerce_map_from_(self, S):
        if isinstance(S, MPolynomialRing_base):
            return True
        return super()._coerce_map_from_(S)

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
            elem = self._from_dict({perm: self.base_ring().one()})
        elif isinstance(x, Permutation):
            perm = x.remove_extra_fixed_points()
            elem = self._from_dict({perm: self.base_ring().one()})
        elif isinstance(x, MPolynomial):
            from sage.interfaces.sympy import sympy_init

            sympy_init()
            sympy_floff = sympify(str(x))
            val = syme.sympify(sympy_floff)
            result = sq.mult_poly(
                {(1, 2): 1},
                val,
                [syme.Symbol(str(g)) for g in self._polynomial_ring.gens()],
                [syme.Symbol(str(g)) for g in self._q_ring.gens()],
            )
            elem = self._from_dict(
                {
                    Permutation(list(k)).remove_extra_fixed_points(): self._q_ring(v)
                    for k, v in result.items()
                }
            )                        
        else:
            elem = None
        
        elem._polynomial_ring = self._polynomial_ring
        elem._q_ring = self._q_ring
        elem._base_varname = self._base_varname
        elem._q_varname = self._q_varname
        return elem

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

def _repr_(self):
    return f"Ring of Quantum Schubert polynomials in {self._base_varname} with {len(self._polynomial_ring.gens())} variables over {self._q_ring.base_ring()}"


FastQuantumSchubertPolynomial = FastQuantumSchubertPolynomial_class
FastQuantumSchubertPolynomialRing_base = FastQuantumSchubertPolynomialRing_xbasis