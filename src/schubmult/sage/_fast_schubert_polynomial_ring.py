from functools import cache

# from sage.misc.parser import Parser
import symengine as syme
from sage.all import *  # noqa: F403
from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.categories.graded_bialgebras_with_basis import GradedBialgebrasWithBasis
from sage.combinat.composition import (
    Composition,
    Compositions,
)
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.permutation import Permutation, Permutations, from_lehmer_code
from sage.misc.cachefunc import cached_method
from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_base
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sympy import sympify

import schubmult.perm_lib as pl
import schubmult.sage._fast_double_schubert_polynomial_ring as bork
import schubmult.schub_lib.double as yz
import schubmult.schub_lib.quantum as sq
import schubmult.schub_lib.quantum_double as qyz
import schubmult.schub_lib.single as py
from schubmult.perm_lib import permtrim

from ._indexing import _coerce_index

# FastQuantumDoubleSchubertPolynomialRing = bork.FastQuantumDoubleSchubertPolynomialRing


def FastSchubertPolynomialRing(
    R: Parent,  # noqa: F405
    num_vars: int,
    base_variable_name: str,
    *,
    code_display: bool = False,
    q_varname: str = "q",
    is_quantum: bool = False,
    indices: tuple[int] = (1,),
):
    """
    Wrapper function to return a double Schubert polynomial Ring

        Calls the _xbasis class to return a (quantum) Schubert
        polynomial ring with the indicated base ring, number of variables,
        variable name, coproduct indices, code_display representation option,
        q-ring variable name, and whether the ring is quantum.

        Example call:

    ```python
    X = FastSchubertPolynomialRing(ZZ, 100, "x")
    X([2, 4, 3, 1]) + X([2, 1, 4, 3])
    ```
    This produces a sum of Schubert polynomials in the "x" variables. These will coerce
    to any polynomial ring with variables with the same names as Schubert polynomials.

    Args:
        R (Parent): The base ring
        num_vars (int): Cardinality of the sets of variables
        base_variable_name (str): Base variable name
        code_display (bool, optional): Whether to display the indices as the Lehmer code. Defaults to False.
        q_varname (str, optional): Variable name of the q-ring. Defaults to "q".
        is_quantum (bool, optional): Whether or not the ring is quantum. Defaults to False.
        indices (tuple[int], optional): Indicies of the variables to split on for the coproduct.

    Returns:
        FastSchubertPolynomialRing_xbasis: Element constructor of the ring

    """
    if is_quantum:
        QR = PolynomialRing(R, num_vars, q_varname)
    else:
        QR = R
    return FastSchubertPolynomialRing_xbasis(
        R,
        num_vars,
        base_variable_name,
        q_varname,
        code_display,
        indices,
        is_quantum,
        QR,
    )


def FastQuantumSchubertPolynomialRing(
    R: Parent,  # noqa: F405
    num_vars: int,
    base_variable_name: str,
    q_varname: str = "q",
    code_display: bool = False,
):
    """
    Quantum Schubert ring generator

    Wraps FastSchubertPolynomialRing(), omitting indices and setting
    is_quantum to True.

    Args:
        R (Parent): The base ring
        num_vars (int): Cardinality of the sets of variables
        base_variable_name (str): Base variable name
        q_varname (str, optional): Variable name of the q-ring. Defaults to "q".
        code_display (bool, optional): Whether to display the indices as the Lehmer code. Defaults to False.


    Returns:
        FastSchubertPolynomialRing_xbasis: Element constructor of the ring

    """
    return FastSchubertPolynomialRing(
        R,
        num_vars,
        base_variable_name,
        q_varname=q_varname,
        code_display=code_display,
        is_quantum=True,
    )


class FastSchubertPolynomial_class(CombinatorialFreeModule.Element):
    @property
    def base_varname(self):
        return self.parent()._base_varname

    @property
    def q_varname(self):
        return self.parent()._q_varname

    @property
    def is_quantum(self):
        return self.parent()._quantum

    @property
    def polynomial_ring(self):
        return self.parent()._polynomial_ring

    def expand(self):
        if self.is_quantum:
            return sum(
                [
                    self.parent()._polynomial_ring(
                        qyz.schubpoly_quantum(
                            pl.Permutation(_coerce_index(k, self.parent()._ascode, False)),
                            self.parent()._polynomial_ring.gens(),
                            [0 for i in range(100)],
                            self.parent()._q_ring.gens(),
                            v,
                        ),
                    )
                    for k, v in self.monomial_coefficients().items()
                ],
            )
        return sum(
            [
                self.parent()._polynomial_ring(
                    yz.schubmult_double(
                        {pl.Permutation([]): v},
                        pl.Permutation(_coerce_index(k, self.parent()._ascode, False)),
                        self.parent()._polynomial_ring.gens(),
                        [0 for i in range(100)],
                    ).get(pl.Permutation([]), 0),
                )
                for k, v in self.monomial_coefficients().items()
            ],
        )

    # PDB STEP WHEN LATEX CALLED
    # def _latex_(self):
    #     rep = "\\mathfrak{S}"
    #     if self.is_quantum:
    #         rep+=f"^{'{'}{self.q_varname}{'}'}"
    #     rep+="_{"+str(()"


@cache
def _single_schub_parser(passed):
    if passed._quantum:
        q_varname = passed._q_varname
    else:
        q_varname = "q"
    QDRing = bork.FastQuantumDoubleSchubertPolynomialRing(
        passed.base_ring(),
        len(passed._polynomial_ring.gens()),
        passed._base_varname,
        coeff_variable_names="y",
        code_display=passed._ascode,
        q_varname=q_varname,
    )
    return QDRing.parser()


class FastSchubertPolynomialRing_xbasis(CombinatorialFreeModule):
    Element = FastSchubertPolynomial_class

    def parser(self):
        if self._parser is None:
            self._parser = _single_schub_parser(self)
        return self._parser

    def __init__(
        self,
        R,
        num_vars,
        base_variable_name,
        q_varname,
        code_display,
        indices,
        quantum,
        QR,
    ):
        self._name = f"{'Quantum ' if quantum else ''}Schubert polynomial ring with X basis"
        self._splitter = indices
        self._repr_option_bracket = False
        self._quantum = quantum

        cat = GradedAlgebrasWithBasis(QR).Commutative() if quantum else GradedBialgebrasWithBasis(R).Commutative()

        index_set = Permutations()
        self._ascode = False

        if code_display:
            index_set = Compositions()
            self._ascode = True

        self._sc_rep = f"{'Q' if quantum else ''}S{base_variable_name}"

        CombinatorialFreeModule.__init__(
            self,
            R if not quantum else QR,
            index_set,
            category=cat,
            prefix=self._sc_rep,
            bracket="(",
        )
        self._q_ring = QR
        self._base_varname = base_variable_name
        self._q_varname = q_varname
        self._polynomial_ring = PolynomialRing(R, num_vars, base_variable_name) if not quantum else PolynomialRing(QR, num_vars, base_variable_name)
        self._populate_coercion_lists_()
        self._parser = None

    def _coerce_map_from_(self, S):
        if isinstance(S, MPolynomialRing_base):
            return True
        if isinstance(S, FastSchubertPolynomialRing_base):
            return True
        if isinstance(S, FastSchubertPolynomialRing_base):
            return True
        if isinstance(S, str):
            return True
        return super()._coerce_map_from_(S)

    @cached_method
    def one_basis(self):
        return _coerce_index([1], False, self._ascode)

    def set_coproduct_indices(self, indices):
        self._splitter = indices

    def _element_constructor_(self, x):
        if isinstance(x, str):
            return self.parser().parse(x)
        if isinstance(x, (list, tuple, Composition, Permutation, pl.Permutation)):
            # checking the input to avoid symmetrica crashing Sage, see trac 12924
            elem = self._from_dict(
                {_coerce_index(x, self._ascode, self._ascode): self.base_ring().one()},
            )
        elif isinstance(x, FastSchubertPolynomial):
            if x.base_varname == self._base_varname and (self._quantum == x.parent()._quantum) and (not self._quantum or x.q_varname == self._q_varname):
                elem = self._from_dict(
                    {_coerce_index(k, x.parent()._ascode, self._ascode): v for k, v in x.monomial_coefficients().items()},
                )
            else:
                return self(x.expand())
        elif isinstance(x, MPolynomial):
            from sage.interfaces.sympy import sympy_init

            sympy_init()
            sympy_floff = sympify(str(x))
            val = syme.sympify(sympy_floff)
            if self._quantum:
                result = sq.mult_poly_q(
                    {pl.Permutation([]): 1},
                    val,
                    [syme.Symbol(str(g)) for g in self._polynomial_ring.gens()],
                    [syme.Symbol(str(g)) for g in self._q_ring.gens()],
                )
            else:
                result = py.mult_poly_py(
                    {pl.Permutation([]): 1},
                    val,
                    [syme.Symbol(str(g)) for g in self._polynomial_ring.gens()],
                )
            elem = self._from_dict(
                {_coerce_index(k, False, self._ascode): self._q_ring(str(v)) if self._quantum else self.base_ring()(str(v)) for k, v in result.items()},
            )
        else:
            raise TypeError(f"Could not convert {x=} to {self}")
        return elem

    def some_elements(self):
        return [
            self.one(),
            self(_coerce_index([1], False, self._ascode)) + 2 * self(_coerce_index([2, 1], False, self._ascode)),
            self(_coerce_index([4, 2, 1, 3], False, self._ascode)) - self(_coerce_index([3, 2, 1], False, self._ascode)),
        ]

    def product_on_basis(self, left, right):
        if self._quantum:
            return sum(
                [
                    self.base_ring()(str(v)) * self(_coerce_index(k, False, self._ascode))
                    for k, v in sq.schubmult_q_fast(
                        {pl.Permutation(_coerce_index(left, self._ascode, False)): 1},
                        pl.Permutation(_coerce_index(right, self._ascode, False)),
                        [syme.sympify(str(g)) for g in self._q_ring.gens()],
                    ).items()
                ],
            )
        return sum(
            [
                self.base_ring()(v) * self(_coerce_index(k, False, self._ascode))
                for k, v in py.schubmult_py(
                    {tuple(_coerce_index(left, self._ascode, False)): 1},
                    tuple(_coerce_index(right, self._ascode, False)),
                ).items()
            ],
        )

    def coproduct_on_basis(self, mperm):
        if self._quantum:
            raise NotImplementedError("Quantum Schubert polynomials do not have a coproduct")
        mperm = _coerce_index(mperm, self._ascode, False)
        indices = self._splitter
        indices = sorted(indices)
        k = len(indices)
        n = len(mperm)
        kcd = [indices[i] - i - 1 for i in range(len(indices))] + [n + 1 - k for i in range(k, n)]
        max_required = max([kcd[i] + i for i in range(len(kcd))])
        kcd2 = kcd + [0 for i in range(len(kcd), max_required)] + [0]
        N = len(kcd)
        kperm = from_lehmer_code(kcd2).inverse()
        coeff_dict = {tuple(kperm): 1}
        coeff_dict = py.schubmult_py(coeff_dict, tuple(mperm))

        inv_kperm = kperm.number_of_inversions()
        inverse_kperm = kperm.inverse()
        total_sum = 0
        for perm, val in coeff_dict.items():
            pperm = Permutation(list(perm))
            downperm = pperm.left_action_product(inverse_kperm)
            if downperm.number_of_inversions() == pperm.number_of_inversions() - inv_kperm:
                flag = True
                for i in range(N):
                    if downperm[i] > N:
                        flag = False
                        break
                if not flag:
                    continue
                firstperm = Permutation(permtrim(list(downperm[0:N])))
                secondperm = Permutation(
                    permtrim([downperm[i] - N for i in range(N, len(downperm))]),
                )
                total_sum += self.base_ring()(val) * self(
                    _coerce_index(firstperm, False, self._ascode),
                ).tensor(self(_coerce_index(secondperm, False, self._ascode)))
        return total_sum

    def _latex_(self, S):
        return f"Boobfucker {S}"

    def _repr_(self):
        return f"Ring of  Schubert polynomials in {self._base_varname} with {len(self._polynomial_ring.gens())} variables over {self._q_ring.base_ring()}"


FastSchubertPolynomial = FastSchubertPolynomial_class
FastSchubertPolynomialRing_base = FastSchubertPolynomialRing_xbasis
