import schubmult.schubmult_yz as yz
from sympy import sympify
import symengine as syme

from ._fast_schubert_polynomial_ring import FastSchubertPolynomialRing_xbasis

from sage.categories.graded_bialgebras_with_basis import GradedBialgebrasWithBasis
from sage.rings.polynomial.flatten import FlatteningMorphism
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.permutation import Permutations, Permutation, from_lehmer_code
from sage.categories.cartesian_product import cartesian_product
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import
from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_base
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

lazy_import("sage.libs.symmetrica", "all", as_="symmetrica")


def FastDoubleSchubertPolynomialRing(
    R, num_vars, varname1, varname2, indices=tuple([1])
):
    """
    Return the FastDoubleSchubert polynomial ring over ``R`` on the X basis.

    This is the basis made of the FastDoubleSchubert polynomials.

    EXAMPLES::

            sage: X = FastDoubleSchubertPolynomialRing(ZZ); X
            Schubert polynomial ring with X basis over Integer Ring
            sage: TestSuite(X).run()
            sage: X(1)
            X[1]
            sage: X([1,2,3])*X([2,1,3])
            X[2, 1]
            sage: X([2,1,3])*X([2,1,3])
            X[3, 1, 2]
            sage: X([2,1,3])+X([3,1,2,4])
            X[2, 1] + X[3, 1, 2]
            sage: a = X([2,1,3])+X([3,1,2,4])
            sage: a^2
            X[3, 1, 2] + 2*X[4, 1, 2, 3] + X[5, 1, 2, 3, 4]
    """

    return FastDoubleSchubertPolynomialRing_xbasis(
        R, num_vars, varname1, varname2, indices
    )


class FastDoubleSchubertPolynomial_class(CombinatorialFreeModule.Element):
    def expand(self):
        """
        EXAMPLES::

                sage: X = FastDoubleSchubertPolynomialRing(ZZ)
                sage: X([2,1,3]).expand()
                x0
                sage: [X(p).expand() for p in Permutations(3)]
                [1, x0 + x1, x0, x0*x1, x0^2, x0^2*x1]

        TESTS:

        Calling .expand() should always return an element of an
        MPolynomialRing::

                sage: X = FastDoubleSchubertPolynomialRing(ZZ)
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
        return sum(
            [
                yz.schubmult(
                    {(1, 2): v},
                    tuple(k[0]),
                    self.parent()._base_polynomial_ring.gens(),
                    self.parent()._coeff_polynomial_rings[k[1]].gens(),
                ).get((1, 2), 0)
                for k, v in self.monomial_coefficients().items()
            ]
        )

    def __eq__(self, other):
        ss = self.parent().one() * self
        oo = self.parent().one() * other
        return ss.monomial_coefficients() == oo.monomial_coefficients()
    
    def __ne__(self, other):
        ss = self.parent().one() * self
        oo = self.parent().one() * other
        return ss.monomial_coefficients() != oo.monomial_coefficients()

    def root_coefficients(self, root_var_name):
        num_vars = len(self.parent()._coeff_polynomial_ring.gens())
        RR = PolynomialRing(
            self.parent().base_ring().base_ring(), num_vars, root_var_name
        )
        r = [sum(RR._first_ngens(j)) for j in range(num_vars)]
        subs_dict = {
            self.parent()._coeff_polynomial_ring.gens()[i]: r[i]
            for i in range(num_vars)
        }
        # res = self
        return self.map_coefficients(lambda foi: RR(foi.subs(subs_dict)))


class FastDoubleSchubertPolynomialRing_xbasis(CombinatorialFreeModule):
    Element = FastDoubleSchubertPolynomial_class

    def __init__(self, R, num_vars, varname1, varname2, indices):
        """
        EXAMPLES::

                sage: X = FastSchubertPolynomialRing(QQ)
                sage: X == loads(dumps(X))
                True
        """
        self._name = "Schubert polynomial ring with X basis"
        self._splitter = indices
        self._repr_option_bracket = False
        self._mixed = False

        if isinstance(varname2, tuple):
            self._mixed = True
            self._varlist = [*varname2]
            self._coeff_polynomial_rings = {
                name: PolynomialRing(R, num_vars, name) for name in self._varlist
            }

            self._coeff_polynomial_ring = R
            for name, CR in self._coeff_polynomial_rings.items():
                self._coeff_polynomial_ring = PolynomialRing(
                    self._coeff_polynomial_ring, num_vars, name
                )
            self._coeff_polynomial_ring = FlatteningMorphism(
                self._coeff_polynomial_ring
            ).codomain()
        else:
            self._varlist = [varname2]
            self._coeff_polynomial_ring = PolynomialRing(R, num_vars, varname2)
            self._coeff_polynomial_rings = {}
            self._coeff_polynomial_rings[varname2] = self._coeff_polynomial_ring

        self._base_polynomial_ring = PolynomialRing(
            self._coeff_polynomial_ring, num_vars, varname1
        )

        self._index_wrapper = cartesian_product([Permutations(), self._varlist])
        cat = GradedBialgebrasWithBasis(self._coeff_polynomial_ring).Commutative()

        CombinatorialFreeModule.__init__(
            self,
            self._coeff_polynomial_ring,
            self._index_wrapper,
            category=cat,
            prefix=f"S({varname1})",
        )
        self._populate_coercion_lists_()

    def set_coproduct_indices(self, indices):
        self._splitter = indices

    

    @cached_method
    def one_basis(self):
        """
        Return the index of the unit of this algebra.

        EXAMPLES::

                sage: X = FastDoubleSchubertPolynomialRing(QQ)
                sage: X.one()  # indirect doctest
                X[1]
        """
        return (Permutation([1]), self._varlist[0])

    def _element_constructor_(self, *x):
        """
        Coerce x into ``self``.

        EXAMPLES::

                sage: X = FastDoubleSchubertPolynomialRing(QQ)
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

                sage: X = FastDoubleSchubertPolynomialRing(QQ)
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
                sage: X = FastDoubleSchubertPolynomialRing(ZZ)
                sage: it = iter(Permutations())
                sage: for _ in range(50):
                ....:     P = next(it)
                ....:     assert X(k(X(P))) == X(P), P
        """
        if len(x) == 1:
            x = x[0]
        elif len(x) > 2:
            raise ValueError("Bad index for element")

        if isinstance(x, list) or isinstance(x, tuple):
            # checking the input to avoid symmetrica crashing Sage, see trac 12924
            if x in self._index_wrapper:
                perm = Permutation(x[0]).remove_extra_fixed_points()
                elem = self._from_dict(
                    {self._index_wrapper((perm, x[1])): self.base_ring().one()}
                )
            elif x in Permutations():
                perm = Permutation(x).remove_extra_fixed_points()
                elem = self._from_dict(
                    {
                        self._index_wrapper(
                            (perm, self._varlist[0])
                        ): self.base_ring().one()
                    }
                )
            elif x[0] in Permutations():
                if x[1] not in self._varlist:
                    raise ValueError(f"{x[1]} is not a valid variable")
                perm = Permutation(x[0]).remove_extra_fixed_points()
                elem = self._from_dict(
                    {self._index_wrapper((perm, x[1])): self.base_ring().one()}
                )
            else:
                raise ValueError("Not perm")
            elem._coeff_polynomial_ring = self._coeff_polynomial_ring
            elem._base_polynomial_ring = self._base_polynomial_ring
            return elem
        elif isinstance(x, Permutation):
            perm = x.remove_extra_fixed_points()
            elem = self._from_dict(
                {self._index_wrapper((perm, self._varlist[0])): self.base_ring().one()}
            )
            elem._coeff_polynomial_ring = self._coeff_polynomial_ring
            elem._base_polynomial_ring = self._base_polynomial_ring
            return elem        
        elif isinstance(x, MPolynomial):
            from sage.interfaces.sympy import sympy_init

            sympy_init()
            sympy_floff = sympify(str(x))
            val = syme.sympify(sympy_floff)
            result = yz.mult_poly(
                {(1, 2): 1},
                val,
                [syme.Symbol(str(g)) for g in self._base_polynomial_ring.gens()],
                [
                    syme.Symbol(str(g))
                    for g in self._coeff_polynomial_rings[self._varlist[0]].gens()
                ],
            )
            elem = self._from_dict(
                {
                    (
                        Permutation(list(k)).remove_extra_fixed_points(),
                        self._varlist[0],
                    ): self._coeff_polynomial_ring(str(v))
                    for k, v in result.items()
                }
            )
            elem._coeff_polynomial_ring = self._coeff_polynomial_ring
            elem._base_polynomial_ring = self._base_polynomial_ring
            return elem
        elif self.has_coerce_map_from(x.parent()):
            try:
                if self._base_polynomial_ring.has_coerce_map_from(x._polynomial_ring):
                    return self(x.expand())
            except Exception:
                raise TypeError
        else:
            raise TypeError


    def some_elements(self):
        """
        Return some elements.

        EXAMPLES::

                sage: X = FastDoubleSchubertPolynomialRing(QQ)
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
                sage: X = FastDoubleSchubertPolynomialRing(QQ)
                sage: X.product_on_basis(p1,p2)
                X[4, 2, 1, 3]
        """
        # return symmetrica.mult_schubert_schubert(left, right)
        # r = [sum(self.base_ring()._first_ngens(j)) for j in range(20)]
        le = tuple(left[0])
        ri = tuple(right[0])
        var_y = [
            self._coeff_polynomial_ring(g)
            for g in self._coeff_polynomial_rings[left[1]].gens()
        ]
        var_z = [
            self._coeff_polynomial_ring(g)
            for g in self._coeff_polynomial_rings[right[1]].gens()
        ]
        result = yz.schubmult(
            {le: 1},
            ri,
            # self._coeff_polynomial_ring.gens(),
            # self._coeff_polynomial_ring.gens(),
            var_y,
            var_z,
        )
        result = {k: v for k, v in result.items() if v != 0}
        return sum(
            [
                self._coeff_polynomial_ring(v)
                * self(Permutation(list(k)), left[1])
                for k, v in result.items()
            ]
        )

    def coproduct_on_basis(self, indm):
        indices = self._splitter
        indices = sorted(indices)
        subs_dict_coprod = {}
        mperm = indm[0]
        mperm = Permutation(list(mperm))
        RR = self._coeff_polynomial_rings[indm[1]]
        RBase = self._coeff_polynomial_rings[self._varlist[0]]
        k = len(indices)
        n = len(mperm)
        kcd = [indices[i] - i - 1 for i in range(len(indices))] + [
            n + 1 - k for i in range(k, n)
        ]
        max_required = max([kcd[i] + i for i in range(len(kcd))])
        kcd2 = kcd + [0 for i in range(len(kcd), max_required)] + [0]
        N = len(kcd)
        kperm = from_lehmer_code(kcd2).inverse()
        # r = [sum(self.base_ring()._first_ngens(j)) for j in range(100)]
        vn = [f"soible_{i}" for i in range(N * 2 + 1)]
        TR = PolynomialRing(self.base_ring(), N * 2 + 1, vn)

        for i in range(N * 2 + 1):
            if i <= N:
                subs_dict_coprod[TR.gens()[i]] = self._coeff_polynomial_ring(
                    RR.gens()[i]
                )
            else:
                subs_dict_coprod[TR.gens()[i]] = self._coeff_polynomial_ring(
                    RBase.gens()[i - N]
                )

        coeff_dict = {tuple(kperm): 1}
        coeff_dict = yz.schubmult(
            coeff_dict,
            tuple(mperm),
            list(TR.gens()),
            RR.gens(),
        )

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
                val = TR(val).subs(subs_dict_coprod)
                total_sum += self._coeff_polynomial_ring(val) * self(
                    (firstperm, indm[1])
                ).tensor(self((secondperm, self._varlist[0])))
        return total_sum

    def _coerce_map_from_(self, S):
        if isinstance(S, MPolynomialRing_base):
            return True
        if isinstance(S, FastSchubertPolynomialRing_xbasis):
            return True
        return super()._coerce_map_from_(S)
        return super()._coerce_map_from_(S)

