from sympy import sympify
import symengine as syme
import schubmult.schubmult_py as py
import schubmult.schubmult_yz as yz
from schubmult.perm_lib import uncode, trimcode, permtrim

from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.permutation import (
    Permutations,
    Permutation,
    from_lehmer_code,
)
from sage.combinat.composition import (
    Compositions,
    Composition,
)
from sage.misc.cachefunc import cached_method
from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_base
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.categories.graded_bialgebras_with_basis import GradedBialgebrasWithBasis

from ._indexing import _coerce_index


def FastSchubertPolynomialRing(
    R, num_vars, varname, indices=tuple([1]), *, code_index=False
):
    return FastSchubertPolynomialRing_xbasis(R, num_vars, varname, indices, code_index)


class FastSchubertPolynomial_class(CombinatorialFreeModule.Element):
    @property
    def base_varname(self):
        return self.parent()._base_varname

    @property
    def polynomial_ring(self):
        return self.parent()._polynomial_ring

    def expand(self):
        return sum(
            [
                yz.schubmult(
                    {(1, 2): v},
                    tuple(_coerce_index(k, self._ascode, False)),
                    self.parent()._polynomial_ring.gens(),
                    [0 for i in range(100)],
                ).get((1, 2), 0)
                for k, v in self.monomial_coefficients().items()
            ]
        )


class FastSchubertPolynomialRing_xbasis(CombinatorialFreeModule):
    Element = FastSchubertPolynomial_class

    def __init__(self, R, num_vars, varname, indices, code_index):
        self._name = "Schubert polynomial ring with X basis"
        self._splitter = indices
        self._repr_option_bracket = False
        cat = GradedBialgebrasWithBasis(R).Commutative()

        index_set = Permutations()
        self._ascode = False

        if code_index:
            index_set = Compositions()
            self._ascode = True

        CombinatorialFreeModule.__init__(
            self,
            R,
            index_set,
            category=cat,
            prefix=f"S{varname}",
        )
        self._polynomial_ring = PolynomialRing(R, num_vars, varname)
        self._base_varname = varname
        self._populate_coercion_lists_()

    def set_coproduct_indices(self, indices):
        self._splitter = indices

    @cached_method
    def one_basis(self):
        return _coerce_index([1], False, self._ascode)

    def _element_constructor_(self, x):
        if isinstance(x, list) or isinstance(x, tuple):
            elem = self._from_dict(
                {_coerce_index(x, self._ascode, self._ascode): self.base_ring().one()}
            )
        elif isinstance(x, FastSchubertPolynomial):
            elem = self._from_dict(
                _coerce_index(
                    x.monomial_coefficients(), x.parent()._ascode, self._ascode
                )
            )
        elif isinstance(x, Permutation):
            elem = self._from_dict(
                {_coerce_index(x, False, self._ascode): self.base_ring().one()}
            )
        elif isinstance(x, Composition):
            elem = self._from_dict(
                {_coerce_index(x, True, self._ascode): self.base_ring().one()}
            )
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
                    _coerce_index(k, False, self._ascode): self.base_ring()(str(v))
                    for k, v in result.items()
                }
            )
        else:
            raise TypeError
        return elem

    def some_elements(self):        
        return [
            self.one(),
            self(_coerce_index([1], False, self._ascode))
            + 2 * self(_coerce_index([2, 1], False, self._ascode)),
            self(_coerce_index([4, 2, 1, 3], False, self._ascode))
            - self(_coerce_index([3, 2, 1], False, self._ascode)),
        ]

    def _coerce_map_from_(self, S):
        if isinstance(S, MPolynomialRing_base):
            return True
        if isinstance(S, FastSchubertPolynomialRing_xbasis):
            return True
        return super()._coerce_map_from_(S)

    def product_on_basis(self, left, right):        
        return sum(
            [
                self.base_ring()(v) * self(_coerce_index(k, False, self._ascode))
                for k, v in py.schubmult(
                    {tuple(_coerce_index(left, self._ascode, False)): 1},
                    tuple(_coerce_index(right, self._ascode, False)),
                ).items()
            ]
        )

    def coproduct_on_basis(self, mperm):
        mperm = _coerce_index(mperm, self._ascode, False)
        indices = self._splitter
        indices = sorted(indices)
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
                total_sum += self.base_ring()(val) * self(
                    _coerce_index(firstperm, False, self._ascode)
                ).tensor(self(_coerce_index(secondperm, False, self._ascode)))
        return total_sum

    def _repr_(self):
        return f"Ring of Schubert polynomials in {self._base_varname} with {len(self._polynomial_ring.gens())} variables over {self.base_ring()} indexed by {'the Lehmer code' if self._ascode else 'permutations'}"


FastSchubertPolynomial = FastSchubertPolynomial_class
FastSchubertPolynomialRing_base = FastSchubertPolynomialRing_xbasis
