from sympy import sympify
import symengine as syme
import schubmult.schubmult_py as py
import schubmult.schubmult_yz as yz

from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.permutation import Permutations, Permutation, from_lehmer_code
from sage.misc.cachefunc import cached_method
# from sage.misc.lazy_import import lazy_import
from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_base
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.categories.graded_bialgebras_with_basis import GradedBialgebrasWithBasis

# lazy_import("sage.libs.symmetrica", "all", as_="symmetrica")


def FastSchubertPolynomialRing(R, num_vars, varname, indices=tuple([1])):
    return FastSchubertPolynomialRing_xbasis(R, num_vars, varname, indices)


class FastSchubertPolynomial_class(CombinatorialFreeModule.Element):
    def expand(self):
        return sum(
            [
                yz.schubmult(
                    {(1, 2): v},
                    tuple(k),
                    self.parent()._polynomial_ring.gens(),
                    [0 for i in range(100)],
                ).get((1, 2), 0)
                for k, v in self.monomial_coefficients().items()
            ]
        )


class FastSchubertPolynomialRing_xbasis(CombinatorialFreeModule):
    Element = FastSchubertPolynomial_class

    def __init__(self, R, num_vars, varname, indices):
        """
        EXAMPLES::

                sage: X = FastSchubertPolynomialRing(QQ)
                sage: X == loads(dumps(X))
                True
        """
        self._name = "Schubert polynomial ring with X basis"
        self._splitter = indices
        self._repr_option_bracket = False
        cat = GradedBialgebrasWithBasis(R)
        CombinatorialFreeModule.__init__(
            self,
            R,
            Permutations(),
            category=cat,
            prefix="S",
            bracket=["_{", "}" + f"({varname})"],
        )
        self._polynomial_ring = PolynomialRing(R, num_vars, varname)
        self._base_varname = varname
        self._populate_coercion_lists_()

    def set_coproduct_indices(self, indices):
        self._splitter = indices

    @cached_method
    def one_basis(self):
        """
        Return the index of the unit of this algebra.

        EXAMPLES::

                sage: X = FastSchubertPolynomialRing(QQ)
                sage: X.one()  # indirect doctest
                X[1]
        """
        return self._indices([1])

    def _element_constructor_(self, x):
        """
        Coerce x into ``self``.

        EXAMPLES::

                sage: X = FastSchubertPolynomialRing(QQ)
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

                sage: X = FastSchubertPolynomialRing(QQ)
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
                sage: X = FastSchubertPolynomialRing(ZZ)
                sage: it = iter(Permutations())
                sage: for _ in range(50):
                ....:     P = next(it)
                ....:     assert X(k(X(P))) == X(P), P
        """
        # print(type(x))
        # elem = None
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
            result = py.mult_poly(
                {(1, 2): 1},
                val,
                [syme.Symbol(str(g)) for g in self._polynomial_ring.gens()],
            )
            elem = self._from_dict(
                {
                    Permutation(list(k)): self.base_ring()(str(v))
                    for k, v in result.items()
                }
            )
        else:
            raise TypeError

        elem._polynomial_ring = self._polynomial_ring
        elem._base_varname = self._base_varname
        return elem

    def some_elements(self):
        """
        Return some elements.

        EXAMPLES::

                sage: X = FastSchubertPolynomialRing(QQ)
                sage: X.some_elements()
                [X[1], X[1] + 2*X[2, 1], -X[3, 2, 1] + X[4, 2, 1, 3]]
        """
        return [
            self.one(),
            self([1]) + 2 * self([2, 1]),
            self([4, 2, 1, 3]) - self([3, 2, 1]),
        ]

    def _coerce_map_from_(self, S):
        if isinstance(S, MPolynomialRing_base):
            return True
        return super()._coerce_map_from_(S)

    def product_on_basis(self, left, right):
        """
        EXAMPLES::

                sage: p1 = Permutation([3,2,1])
                sage: p2 = Permutation([2,1,3])
                sage: X = FastSchubertPolynomialRing(QQ)
                sage: X.product_on_basis(p1,p2)
                X[4, 2, 1, 3]
        """
        # return symmetrica.mult_schubert_schubert(left, right)
        return sum(
            [
                self.base_ring()(v) * self(Permutation(list(k)))
                for k, v in py.schubmult({tuple(left): 1}, tuple(right)).items()
            ]
        )

    def coproduct_on_basis(self, mperm):
        """
        Coproduct on a single Schubert polynomial
        Depends on the indices.
        """
        indices = self._splitter
        indices = sorted(indices)
        mperm = Permutation(list(mperm))
        k = len(indices)
        n = len(mperm)
        kcd = [indices[i] - i - 1 for i in range(len(indices))] + [
            n + 1 - k for i in range(k, n)
        ]
        max_required = max([kcd[i] + i for i in range(len(kcd))])
        kcd2 = kcd + [0 for i in range(len(kcd), max_required)] + [0]
        N = len(kcd)
        kperm = from_lehmer_code(kcd2).inverse()
        coeff_dict = {tuple(kperm): 1}
        coeff_dict = py.schubmult(coeff_dict, tuple(mperm))

        inv_kperm = kperm.number_of_inversions()
        inverse_kperm = kperm.inverse()
        total_sum = 0
        for perm, val in coeff_dict.items():
            pperm = Permutation(list(perm))
            downperm = pperm.left_action_product(inverse_kperm)
            if (
                downperm.number_of_inversions()
                == pperm.number_of_inversions() - inv_kperm
            ):
                flag = True
                for i in range(N):
                    if downperm[i] > N:
                        flag = False
                        break
                if not flag:
                    continue
                firstperm = Permutation(list(downperm[0:N]))
                secondperm = Permutation(
                    [downperm[i] - N for i in range(N, len(downperm))]
                )
                total_sum += self.base_ring()(val) * self(firstperm).tensor(
                    self(secondperm)
                )
        return total_sum
    
    def _repr_(self):
        return f"Ring of Schubert polynomials in {self._base_varname} with {len(self._polynomial_ring.gens())} variables"



FastSchubertPolynomial = FastSchubertPolynomial_class
FastSchubertPolynomialRing_base = FastSchubertPolynomialRing_xbasis