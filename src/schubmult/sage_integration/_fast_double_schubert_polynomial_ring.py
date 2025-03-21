import schubmult.schubmult_q_yz as qyz
import schubmult.schubmult_yz as yz
from sympy import sympify
import symengine as syme

# from ._fast_quantum_schubert_polynomial_ring import (
#     FastSchubertPolynomial,
#     FastSchubertPolynomialRing_base,
# )

from ._indexing import _coerce_index
from sage.combinat.composition import (
    Compositions,
    Composition,
)

from ._fast_schubert_polynomial_ring import (
    FastSchubertPolynomialRing_base,
    FastSchubertPolynomial,
)

from ._fast_double_schubert_polynomial_ring import (
    FastDoubleSchubertPolynomialRing_base,
    FastDoubleSchubertPolynomial,
)

from sage.categories.graded_bialgebras_with_basis import GradedBialgebrasWithBasis
from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.categories.cartesian_product import cartesian_product
from sage.combinat.permutation import Permutations, Permutation, from_lehmer_code
from sage.misc.cachefunc import cached_method
from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_base
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.flatten import FlatteningMorphism


def FastDoubleSchubertPolynomialRing(
    R,
    num_vars,
    varname1,
    varname2,
    *,
    indices=tuple([1]),
    code_index=False,
    q_varname="q",
    quantum=False,
):
    QR = None
    if quantum:
        QR = PolynomialRing(R, num_vars, q_varname)
    return FastDoubleSchubertPolynomialRing_xbasis(
        R,
        num_vars,
        varname1,
        varname2,
        q_varname,
        quantum,
        code_index,
        indices,
        quantum,
        QR,
    )


def FastDoubleQuantumSchubertPolynomialRing(
    R,
    num_vars,
    varname1,
    varname2,
    *,
    code_index=False,
    q_varname="q",
):
    return FastDoubleSchubertPolynomialRing(
        R,
        num_vars,
        varname1,
        varname2,
        code_index=code_index,
        indices=tuple([1]),
        quantum=True,
        q_varname=q_varname,
    )


class FastDoubleSchubertPolynomial_class(CombinatorialFreeModule.Element):
    @property
    def base_varname(self):
        return self.parent()._base_varname

    @property
    def q_varname(self):
        if not self.parent().is_quantum:
            return None
        else:
            return self.parent()._q_varname

    @property
    def base_polynomial_ring(self):
        return self.parent()._base_polynomial_ring

    @property
    def coeff_polynomial_ring(self):
        return self.parent()._coeff_polynomial_ring

    @property
    def q_ring(self):
        if not self.parent()._quantum:
            return None
        else:
            return self.parent()._q_ring

    def expand(self):
        if self.parent()._quantum:
            return sum(
                [
                    qqyz.schubpoly_quantum(
                        tuple(_coerce_index(k[0], self.parent()._ascode, False)),
                        self.parent()._base_polynomial_ring.gens(),
                        self.parent()._coeff_polynomial_rings[k[1]].gens(),
                        self.parent()._q_ring.gens(),
                        v,
                    )
                    for k, v in self.monomial_coefficients().items()
                ]
            )
        else:
            return sum(
                [
                    yz.schubmult(
                        {(1, 2): v},
                        tuple(_coerce_index(k[0], self.parent()._ascode, False)),
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

    def __init__(
        self,
        R,
        num_vars,
        varname1,
        varname2,
        q_varname,
        code_index,
        indices,
        quantum,
        QR,
    ):
        self._name = (
            f"{'Quantum double' if quantum else 'Double'} Schubert polynomial ring"
        )
        self._splitter = indices
        self._repr_option_bracket = False
        self._mixed = False
        self._q_ring = QR
        self._quantum = quantum
        self._base_varname = varname1
        self._q_varname = q_varname

        if isinstance(varname2, tuple):
            self._mixed = True
            self._varlist = [*varname2]
            self._coeff_polynomial_rings = {
                name: PolynomialRing(R if not quantum else QR, num_vars, name)
                for name in self._varlist
            }

            self._coeff_polynomial_ring = R if not quantum else QR
            for name, CR in self._coeff_polynomial_rings.items():
                self._coeff_polynomial_ring = PolynomialRing(
                    self._coeff_polynomial_ring, num_vars, name
                )
            self._coeff_polynomial_ring = FlatteningMorphism(
                self._coeff_polynomial_ring
            ).codomain()
        else:
            self._varlist = [varname2]
            self._coeff_polynomial_ring = PolynomialRing(
                R if not quantum else QR, num_vars, varname2
            )
            self._coeff_polynomial_rings = {}
            self._coeff_polynomial_rings[varname2] = self._coeff_polynomial_ring

        index_set = Permutations()
        self._ascode = False

        if code_index:
            index_set = Compositions()
            self._ascode = True

        self._base_polynomial_ring = PolynomialRing(
            self._coeff_polynomial_ring, num_vars, varname1
        )

        self._index_wrapper = cartesian_product([index_set, self._varlist])
        cat = (
            GradedAlgebrasWithBasis(self._coeff_polynomial_ring).Commutative()
            if quantum
            else GradedBialgebrasWithBasis(self._coeff_polynomial_ring).Commutative()
        )

        CombinatorialFreeModule.__init__(
            self,
            self._coeff_polynomial_ring,
            self._index_wrapper,
            category=cat,
            prefix=f"{'Q' if quantum else ''}S{varname1}",
        )
        self._populate_coercion_lists_()

    @cached_method
    def one_basis(self):
        return (_coerce_index([1], False, self._ascode), self._varlist[0])

    @property
    def is_quantum(self):
        return self._quantum

    def _element_constructor_(self, *x):
        if len(x) == 1:
            x = x[0]
        elif len(x) > 2:
            raise ValueError("Bad index for element")

        if (
            isinstance(x, list)
            or isinstance(x, tuple)
            or isinstance(x, Permutation)
            or isinstance(x, Composition)
        ):
            # checking the input to avoid symmetrica crashing Sage, see trac 12924
            if x in self._index_wrapper:
                elem = self._from_dict(
                    {self._index_wrapper((x[0], x[1])): self.base_ring().one()}
                )
            else:
                elem = self._from_dict(
                    {
                        self._index_wrapper(
                            (
                                _coerce_index(x, self._ascode, self._ascode),
                                self._varlist[0],
                            )
                        ): self.base_ring().one()
                    }
                )
        elif isinstance(x, MPolynomial):
            from sage.interfaces.sympy import sympy_init

            sympy_init()
            sympy_floff = sympify(str(x))
            val = syme.sympify(sympy_floff)
            if self._quantum:
                result = qyz.mult_poly(
                    {(1, 2): 1},
                    val,
                    [syme.Symbol(str(g)) for g in self._base_polynomial_ring.gens()],
                    [
                        syme.Symbol(str(g))
                        for g in self._coeff_polynomial_rings[self._varlist[0]].gens()
                    ],
                    [syme.Symbol(str(g)) for g in self._q_ring.gens()],
                )
            else:
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
                        _coerce_index(k, False, self._ascode),
                        self._varlist[0],
                    ): self._coeff_polynomial_ring(str(v))
                    for k, v in result.items()
                }
            )
        elif isinstance(x, FastSchubertPolynomial) or isinstance(
            x, FastDoubleSchubertPolynomial
        ):
            return self(x.expand())
        elif isinstance(x, FastDoubleSchubertPolynomial):
            if (
                x.base_varname == self._base_varname
                and (self._quantum == x.parent()._quantum)
                and (not self._quantum or x.q_varname == self._q_varname)
            ):
                elem = self._from_dict(
                    {
                        (_coerce_index(k[0], x.parent()._ascode, self._ascode), k[1]): v
                        for k, v in x.monomial_coefficients().items()
                    }
                )
            else:
                return self(x.expand())
        elif isinstance(x, FastSchubertPolynomial):
            if (
                x.base_varname == self._base_varname
                and x.q_varname == self._q_varname
                and (self._quantum == x.parent()._quantum)
                and (not self._quantum or x.q_varname == self._q_varname)
            ):
                elem_dict = {}

                for k, v in x.monomial_coefficients().items():
                    if self._quantum:
                        res = qyz.schubmult_db(
                            {(1, 2): self._coeff_polynomial_ring(v)},
                            tuple(_coerce_index(k, x.parent()._ascode, False)),
                            self._coeff_polynomial_rings[self._varlist[0]].gens(),
                            [
                                0
                                for i in range(
                                    len(
                                        self._coeff_polynomial_rings[
                                            self._varlist[0]
                                        ].gens()
                                    )
                                )
                            ],
                            self._q_ring.gens(),
                        )
                    else:
                        res = yz.schubmult(
                            {(1, 2): self._coeff_polynomial_ring(v)},
                            tuple(_coerce_index(k, x.parent()._ascode, False)),
                            self._coeff_polynomial_rings[self._varlist[0]].gens(),
                            [
                                0
                                for i in range(
                                    len(
                                        self._coeff_polynomial_rings[
                                            self._varlist[0]
                                        ].gens()
                                    )
                                )
                            ],
                        )
                    for k0, c0 in res.items():
                        elem_dict[
                            (_coerce_index(k0, False, self._ascode), self._varlist[0])
                        ] = elem_dict.get(
                            (_coerce_index(k0, False, self._ascode), self._varlist[0]),
                            self._coeff_polynomial_ring.zero(),
                        ) + self._coeff_polynomial_ring(c0)
                elem = self._from_dict(elem_dict)
            else:
                elem = self(x.expand())
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

    def product_on_basis(self, left, right):
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
        if self._quantum:
            result = qyz.schubmult_db(
                {tuple(_coerce_index(le, self._ascode, False)): 1},
                tuple(_coerce_index(ri, self._ascode, False)),
                var_y,
                var_z,
                self._q_ring.gens(),
            )
        else:
            result = yz.schubmult(
                {tuple(_coerce_index(le, self._ascode, False)): 1},
                tuple(_coerce_index(ri, self._ascode, False)),
                var_y,
                var_z,
            )
        result = {k: v for k, v in result.items() if v != 0}
        return sum(
            [
                self._coeff_polynomial_ring(v)
                * self((_coerce_index(k, False, self._ascode), left[1]))
                for k, v in result.items()
            ]
        )

    def _coerce_map_from_(self, S):
        if isinstance(S, MPolynomialRing_base):
            return True
        if isinstance(S, FastSchubertPolynomialRing_base):
            return True
        if isinstance(S, FastDoubleSchubertPolynomialRing_base):
            return True
        return super().has_coerce_map_from(S)

    def coproduct_on_basis(self, indm):
        if self._quantum:
            raise NotImplementedError("Quantum double Schubert polynomials have no coproduct")
        indices = self._splitter
        indices = sorted(indices)
        subs_dict_coprod = {}
        mperm = indm[0]
        mperm = _coerce_index(mperm, self._ascode, False)
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
                    (_coerce_index(firstperm, False, self._ascode), indm[1])
                ).tensor(
                    self(
                        (
                            _coerce_index(secondperm, False, self._ascode),
                            self._varlist[0],
                        )
                    )
                )
        return total_sum

    def _repr_(self):
        return (
            f"Ring of {'quantum' if self._quantum else ''} double Schubert polynomials in {self._base_varname}{',' + self._q_varname if self._quantum else ''} with {len(self._base_polynomial_ring.gens())} variables with"
            f" coefficient variables {','.join(self._varlist)} over the ring {self._coeff_polynomial_ring.base_ring()} indexed by {'the Lehmer code' if self._ascode else 'permutations'}"
        )


FastDoubleSchubertPolynomial = FastDoubleSchubertPolynomial_class
FastDoubleSchubertPolynomialRing_base = FastDoubleSchubertPolynomialRing_xbasis
