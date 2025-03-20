import schubmult.schubmult_q_yz as yz
from sympy import sympify
import symengine as syme

from ._fast_quantum_schubert_polynomial_ring import (
    FastQuantumSchubertPolynomial,
    FastQuantumSchubertPolynomialRing_base
)

from ._fast_schubert_polynomial_ring import FastSchubertPolynomialRing_base, FastSchubertPolynomial

from ._fast_double_schubert_polynomial_ring import FastDoubleSchubertPolynomialRing_base, FastDoubleSchubertPolynomial


from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.categories.cartesian_product import cartesian_product
from sage.combinat.permutation import Permutations, Permutation
from sage.misc.cachefunc import cached_method
from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_base
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.flatten import FlatteningMorphism


def FastQuantumDoubleSchubertPolynomialRing(
    R, num_vars, varname1, varname2, q_varname="q_"
):
    """
    Return the FastQuantumDoubleSchubert polynomial ring over ``R`` on the X basis.

    This is the basis made of the FastQuantumDoubleSchubert polynomials.

    EXAMPLES::

            sage: X = FastQuantumDoubleSchubertPolynomialRing(ZZ); X
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
    QR = PolynomialRing(R, num_vars, q_varname)
    return FastQuantumDoubleSchubertPolynomialRing_xbasis(
        QR, num_vars, varname1, varname2, q_varname
    )


class FastQuantumDoubleSchubertPolynomial_class(CombinatorialFreeModule.Element):

    @property
    def base_varname(self):
        return self.parent()._base_varname
    
    @property
    def q_varname(self):
        return self.parent()._q_varname

    @property
    def base_polynomial_ring(self):
        return self.parent()._base_polynomial_ring
    
    
    @property
    def coeff_polynomial_ring(self):
        return self.parent()._coeff_polynomial_ring
    
    
    @property
    def q_ring(self):
        return self.parent()._q_ring


    def expand(self):
        return sum(
            [
                yz.schubpoly_quantum(
                    tuple(k[0]),
                    self.parent()._base_polynomial_ring.gens(),
                    self.parent()._coeff_polynomial_rings[k[1]].gens(),
                    self.parent()._q_ring.gens(),
                    v,
                )
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


class FastQuantumDoubleSchubertPolynomialRing_xbasis(CombinatorialFreeModule):
    Element = FastQuantumDoubleSchubertPolynomial_class

    def __init__(self, R, num_vars, varname1, varname2, q_varname):
        """
        EXAMPLES::

                sage: X = FastQuantumDoubleSchubertPolynomialRing(QQ)
                sage: X == loads(dumps(X))
                True
        """
        self._name = "QuantumDouble Schubert polynomial ring with X basis"
        self._repr_option_bracket = False
        self._mixed = False
        self._q_ring = R
        self._base_varname = varname1
        self._q_varname = q_varname

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
        cat = GradedAlgebrasWithBasis(self._coeff_polynomial_ring).Commutative()

        CombinatorialFreeModule.__init__(
            self,
            self._coeff_polynomial_ring,
            self._index_wrapper,
            category=cat,
            prefix=f"S^{q_varname}({varname1})",
        )
        self._populate_coercion_lists_()


    @cached_method
    def one_basis(self):
        """
        Return the index of the unit of this algebra.

        EXAMPLES::

                sage: X = FastQuantumDoubleSchubertPolynomialRing(QQ)
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
                raise ValueError            
        elif isinstance(x, Permutation):
            perm = x.remove_extra_fixed_points()
            elem = self._from_dict(
                {self._index_wrapper((perm, self._varlist[0])): self.base_ring().one()}
            )            
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
                [syme.Symbol(str(g)) for g in self._q_ring.gens()],
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
        elif isinstance(x, FastSchubertPolynomial) or isinstance(x, FastDoubleSchubertPolynomial):
            return self(x.expand())
        elif isinstance(x, FastQuantumSchubertPolynomial):
            import inspect
            if x.base_varname == self._base_varname and x.q_varname == self._q_varname:
                elem_dict = {}
                for k, v in x.monomial_coefficients().items():
                    res = yz.schubmult_db(
                                {(1, 2): self._coeff_polynomial_ring(v)},
                                tuple(k),
                                self._coeff_polynomial_rings[self._varlist[0]].gens(),
                                [
                                    0
                                    for i in range(
                                        len(self._coeff_polynomial_ring.gens())
                                    )
                                ],
                                self._q_ring.gens()
                            )
                    for k0, c0 in res.items():
                        elem_dict[(k0, self._varlist[0])] = elem_dict.get((k0, self._varlist[0]), self._coeff_polynomial_ring.zero()) + self._coeff_polynomial_ring(c0)                        
                elem = self._from_dict(elem_dict)
            else:
                elem = self(x.expand())
        else:
            elem = None         

        
        return elem

    def some_elements(self):
        """
        Return some elements.

        EXAMPLES::

                sage: X = FastQuantumDoubleSchubertPolynomialRing(QQ)
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
                sage: X = FastQuantumDoubleSchubertPolynomialRing(QQ)
                sage: X.product_on_basis(p1,p2)
                X[4, 2, 1, 3]
        """
        # return symmetrica.mult_schubert_schubert(left, right)
        # r = [sum(self.base_ring()._first_ngens(j)) for j in range(100)]
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
        result = yz.schubmult_db(
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
                * self((Permutation(list(k)), left[1]))
                for k, v in result.items()
            ]
        )

    def _coerce_map_from_(self, S):
        if isinstance(S, MPolynomialRing_base):
            return True
        if isinstance(S, FastQuantumSchubertPolynomialRing_base):
            return True
        if isinstance(S, FastSchubertPolynomialRing_base):
            return True
        if isinstance(S, FastDoubleSchubertPolynomialRing_base):
            return True
        return super().has_coerce_map_from(S)
    
    def _repr_(self):
        return (
            f"Ring of quantum double Schubert polynomials in {self._base_varname},{self._q_varname} with {len(self._base_polynomial_ring.gens())} variables with"
            f" coefficient variables {','.join(self._varlist)} over the ring {self._coeff_polynomial_ring.base_ring()}"
        )
FastQuantumDoubleSchubertPolynomial = FastQuantumDoubleSchubertPolynomial_class
FastQuantumDoubleSchubertPolynomialRing_base = FastQuantumDoubleSchubertPolynomialRing_xbasis