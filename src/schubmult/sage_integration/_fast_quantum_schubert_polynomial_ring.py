import schubmult.schubmult_q as sq
import schubmult.schubmult_q_yz as qyz
from ._indexing import _coerce_index

from sympy import sympify
import symengine as syme

from ._fast_schubert_polynomial_ring import (
    FastSchubertPolynomial,
    FastSchubertPolynomialRing_base,
)

from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.combinat.free_module import CombinatorialFreeModule

from sage.combinat.permutation import Permutations, Permutation
from sage.combinat.composition import (
    Compositions,
    Composition,
)
from sage.misc.cachefunc import cached_method
from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_base
from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing


def FastQuantumSchubertPolynomialRing(
    R, num_vars, varname, q_varname="q", *, code_index=False
):
    QR = PolynomialRing(R, num_vars, q_varname)
    return FastQuantumSchubertPolynomialRing_xbasis(
        QR, num_vars, varname, q_varname, code_index
    )


class FastQuantumSchubertPolynomial_class(CombinatorialFreeModule.Element):
    @property
    def base_varname(self):
        return self.parent()._base_varname

    @property
    def q_varname(self):
        return self.parent()._q_varname

    @property
    def polynomial_ring(self):
        return self.parent()._polynomial_ring

    def expand(self):
        return sum(
            [
                self.parent()._polynomial_ring(
                    qyz.schubpoly_quantum(
                        tuple(_coerce_index(k, self.parent()._ascode, False)),
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

    def __init__(self, R, num_vars, varname, q_varname, code_index):
        self._name = "QuantumSchubert polynomial ring with X basis"
        self._repr_option_bracket = False
        cat = GradedAlgebrasWithBasis(
            R
        ).Commutative()  # CoalgebrasWithBasis(R).Graded()

        index_set = Permutations()
        self._ascode = False

        if code_index:
            index_set = Compositions()
            self._ascode = True

        CombinatorialFreeModule.__init__(
            self, R, index_set, category=cat, prefix=f"QS{varname}"
        )
        self._q_ring = R
        self._base_varname = varname
        self._q_varname = q_varname
        self._polynomial_ring = PolynomialRing(R, num_vars, varname)
        self._populate_coercion_lists_()

    def _coerce_map_from_(self, S):
        if isinstance(S, MPolynomialRing_base):
            return True
        if isinstance(S, FastSchubertPolynomialRing_base):
            return True
        if isinstance(S, FastQuantumSchubertPolynomialRing_base):
            return True
        return super()._coerce_map_from_(S)

    @cached_method
    def one_basis(self):
        return _coerce_index([1], False, self._ascode)

    def _element_constructor_(self, x):
        if (
            isinstance(x, list)
            or isinstance(x, tuple)
            or isinstance(x, Composition)
            or isinstance(x, Permutation)
        ):
            # checking the input to avoid symmetrica crashing Sage, see trac 12924
            elem = self._from_dict(
                {_coerce_index(x, self._ascode, self._ascode): self.base_ring().one()}
            )
        elif isinstance(x, FastQuantumSchubertPolynomial):
            elem = self._from_dict(
                _coerce_index(
                    x.monomial_coefficients(), x.parent()._ascode, self._ascode
                )
            )
        elif isinstance(x, MPolynomial) or isinstance(x, FastSchubertPolynomial):
            if isinstance(x, FastSchubertPolynomial):
                x = x.expand()
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
                    _coerce_index(k, False, self._ascode): self._q_ring(v)
                    for k, v in result.items()
                }
            )
        else:
            raise TypeError(f"Could not convert {x=} to {self}")
        return elem

    def some_elements(self):
        return [
            self.one(),
            self(_coerce_index([1], False, self._ascode))
            + 2 * self(_coerce_index([2, 1], False, self._ascode)),
            self(_coerce_index([4, 2, 1, 3], False, self._ascode))
            - self(_coerce_index([3, 2, 1], False, self._ascode)),
        ]

    def product_on_basis(self, left, right):
        return sum(
            [
                self.base_ring()(v) * self(_coerce_index(k, False, self._ascode))
                for k, v in sq.schubmult_db(
                    {
                        tuple(
                            _coerce_index(left, self._ascode, False)
                        ): self.base_ring()(1)
                    },
                    tuple(_coerce_index(right, self._ascode, False)),
                    list(self.base_ring().gens()),
                ).items()
            ]
        )


def _repr_(self):
    return f"Ring of Quantum Schubert polynomials in {self._base_varname} with {len(self._polynomial_ring.gens())} variables over {self._q_ring.base_ring()}"


FastQuantumSchubertPolynomial = FastQuantumSchubertPolynomial_class
FastQuantumSchubertPolynomialRing_base = FastQuantumSchubertPolynomialRing_xbasis
