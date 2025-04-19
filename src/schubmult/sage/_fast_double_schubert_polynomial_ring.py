from functools import cache

import symengine as syme
from sage.all import *  # noqa: F403
from sage.categories.cartesian_product import cartesian_product
from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.categories.graded_bialgebras_with_basis import GradedBialgebrasWithBasis
from sage.combinat.composition import (
    Composition,
    Compositions,
)
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.permutation import Permutation, Permutations
from sage.misc.cachefunc import cached_method
from sage.misc.parser import Parser
from sage.rings.polynomial.flatten import FlatteningMorphism
from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_base
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sympy import sympify

import schubmult.perm_lib as pl

# from . import (
#     FastSchubertPolynomialRing_base,
#     FastSchubertPolynomialRing,
#     FastSchubertPolynomial,
# )as
import schubmult.sage._fast_schubert_polynomial_ring as bork
import schubmult.schub_lib.double as yz
import schubmult.schub_lib.quantum_double as qyz
from schubmult.perm_lib import permtrim

from ._indexing import _coerce_index


def FastDoubleSchubertPolynomialRing(
    R: Parent,  # noqa: F405
    num_vars: int,
    base_variable_name: str,
    coeff_variable_names: str | tuple[str],
    *,
    indices: tuple[int] = (1,),
    code_display: bool = False,
    q_varname: str = "q",
    is_quantum: bool = False,
):
    """
    Wrapper function to return a double Schubert polynomial Ring

        Calls the _xbasis class to return a double or quantum double Schubert
        polynomial ring with the indicated base ring, number of variables,
        variable names (base variable, and then one or more sets of coefficient)
        variables, coproduct indices, code_display representation option, q-ring
        variable name, and whether the ring is quantum.

        Example call:

    ```python
    X = FastDoubleSchubertPolynomialRing(ZZ, 100, "x", ("y", "z"))
    X([2, 4, 3, 1]) + X([2, 1, 4, 3], "z")
    ```

    Args:
            R (sage ring): The base ring
            num_vars (int): Cardinality of the sets of variables
            base_variable_name (str): Base variable name
            coeff_variable_names (str | tuple[str]): Coefficient variable name(s)
            indices (tuple[int], optional): Indicies of the variables to split on for the coproduct.
            code_display (bool, optional): Whether to display the indices as the Lehmer code. Defaults to False.
            q_varname (str, optional): Variable name of the q-ring. Defaults to "q".
            is_quantum (bool, optional): Whether or not the ring is quantum. Defaults to False.

    Returns:
            FastDoubleSchubertPolynomialRing_xbasis: Basis element generator of the ring

    """
    QR = None
    if is_quantum:
        QR = PolynomialRing(R, num_vars, q_varname)
    return FastDoubleSchubertPolynomialRing_xbasis(
        R,
        num_vars,
        base_variable_name,
        coeff_variable_names,
        q_varname,
        code_display,
        indices,
        is_quantum,
        QR,
    )


def FastQuantumDoubleSchubertPolynomialRing(
    R: Parent,  # noqa: F405
    num_vars: int,
    base_variable_name: str,
    coeff_variable_names: str | tuple[str],
    *,
    code_display=False,
    q_varname="q",
):
    """
    Quantum double Schubert ring generator

    Wraps FastDoubleSchubertPolynomialRing(), omitting indices and setting
    is_quantum to True.

    Args:
        R (sage ring): The base ring
        num_vars (int): Cardinality of the sets of variables
        base_variable_name (str): Base variable name
        coeff_variable_names (str | tuple[str]): Coefficient variable name(s)
        code_display (bool, optional): Whether to display the indices as the Lehmer code. Defaults to False.
        q_varname (str, optional): Variable name of the q-ring. Defaults to "q".

    Returns:
        FastDoubleSchubertPolynomialRing_xbasis: Basis element generator of the quantum ring

    """
    return FastDoubleSchubertPolynomialRing(
        R,
        num_vars,
        base_variable_name,
        coeff_variable_names,
        code_display=code_display,
        indices=(1,),
        is_quantum=True,
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
        return self.parent()._q_ring

    # can speed this up
    def expand(self):
        if self.parent()._quantum:
            return sum(
                [
                    qyz.schubpoly_quantum(
                        pl.Permutation(_coerce_index(k[0], self.parent()._ascode, False)),
                        self.parent()._base_polynomial_ring.gens(),
                        self.parent()._coeff_polynomial_rings[k[1]].gens(),
                        self.parent()._q_ring.gens(),
                        v,
                    )
                    for k, v in self.monomial_coefficients().items()
                ],
            )
        return sum(
            [
                yz.schubmult_double(
                    {pl.Permutation([]): v},
                    pl.Permutation(_coerce_index(k[0], self.parent()._ascode, False)),
                    self.parent()._base_polynomial_ring.gens(),
                    self.parent()._coeff_polynomial_rings[k[1]].gens(),
                ).get(pl.Permutation([]), 0)
                for k, v in self.monomial_coefficients().items()
            ],
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
        RR = PolynomialRing(self.parent().base_ring().base_ring(), num_vars, root_var_name)
        r = [sum(RR._first_ngens(j)) for j in range(num_vars)]
        subs_dict = {self.parent()._coeff_polynomial_ring.gens()[i]: r[i] for i in range(num_vars)}
        # res = self
        return self.map_coefficients(lambda foi: RR(foi.subs(subs_dict)))


@cache
def _double_schub_parser(passed):
    fdict = {}
    vardict = {}
    fdict[passed._sc_rep] = passed._element_constructor_
    if passed._quantum:
        QSRing = bork.FastSchubertPolynomialRing(
            passed._base_polynomial_ring.base_ring(),
            len(passed._base_polynomial_ring.gens()),
            passed._base_varname,
            is_quantum=True,
            code_display=passed._ascode,
            q_varname=passed._q_varname,
        )
        fdict[QSRing._sc_rep] = QSRing._element_constructor_
    SRing = bork.FastSchubertPolynomialRing(
        passed._base_polynomial_ring.base_ring().base_ring(),
        len(passed._base_polynomial_ring.gens()),
        passed._base_varname,
        is_quantum=False,
        code_display=passed._ascode,
    )
    DRing = FastDoubleSchubertPolynomialRing(
        passed._base_polynomial_ring.base_ring().base_ring(),
        len(passed._base_polynomial_ring.gens()),
        passed._base_varname,
        coeff_variable_names=tuple(passed._varlist),
    )
    fdict[DRing._sc_rep] = DRing._element_constructor_
    fdict[SRing._sc_rep] = SRing._element_constructor_
    for g in passed._base_polynomial_ring.gens():
        vardict[str(g)] = g
    return Parser(make_function=fdict, make_var=vardict)


class FastDoubleSchubertPolynomialRing_xbasis(CombinatorialFreeModule):
    Element = FastDoubleSchubertPolynomial_class

    # def inject_variables:

    def parser(self):
        if self._parser is None:
            self._parser = _double_schub_parser(self)
        return self._parser

    def __init__(
        self,
        R,
        num_vars,
        base_variable_name,
        coeff_variable_names,
        q_varname,
        code_display,
        indices,
        quantum,
        QR,
    ):
        self._name = f"{'Quantum double' if quantum else 'Double'} Schubert polynomial ring"
        self._splitter = indices
        self._repr_option_bracket = False
        self._mixed = False
        self._q_ring = QR
        self._quantum = quantum
        self._base_varname = base_variable_name
        self._q_varname = q_varname

        if isinstance(coeff_variable_names, tuple):
            self._mixed = True
            self._varlist = [*coeff_variable_names]
            self._coeff_polynomial_rings = {name: PolynomialRing(R if not quantum else QR, num_vars, name) for name in self._varlist}

            self._coeff_polynomial_ring = R if not quantum else QR
            for name, CR in self._coeff_polynomial_rings.items():
                self._coeff_polynomial_ring = PolynomialRing(
                    self._coeff_polynomial_ring,
                    num_vars,
                    name,
                )
            self._coeff_polynomial_ring = FlatteningMorphism(self._coeff_polynomial_ring).codomain()
        else:
            self._varlist = [coeff_variable_names]
            self._coeff_polynomial_ring = PolynomialRing(
                R if not quantum else QR,
                num_vars,
                coeff_variable_names,
            )
            self._coeff_polynomial_rings = {}
            self._coeff_polynomial_rings[coeff_variable_names] = self._coeff_polynomial_ring

        index_set = Permutations()
        self._ascode = False

        if code_display:
            index_set = Compositions()
            self._ascode = True

        self._base_polynomial_ring = PolynomialRing(
            self._coeff_polynomial_ring,
            num_vars,
            base_variable_name,
        )

        self._index_wrapper = cartesian_product([index_set, self._varlist])
        cat = GradedAlgebrasWithBasis(self._coeff_polynomial_ring).Commutative() if quantum else GradedBialgebrasWithBasis(self._coeff_polynomial_ring).Commutative()

        self._sc_rep = f"{'Q' if quantum else ''}DS{base_variable_name}"

        #         [('bracket', None),
        #  ('iterate_key', False),
        #  ('latex_bracket', False), ('latex_names', None),
        #  ('latex_prefix', None), ('latex_scalar_mult', None),
        #  ('names', None), ('prefix', 'x'),
        #  ('scalar_mult', '*'),
        #  ('sorting_key', <function ...<lambda> at ...>),
        #  ('sorting_reverse', False), ('string_quotes', True),
        #  ('tensor_symbol', None)]

        CombinatorialFreeModule.__init__(
            self,
            self._coeff_polynomial_ring,
            self._index_wrapper,
            category=cat,
            prefix=self._sc_rep,
        )
        # print_options
        self._populate_coercion_lists_()
        self._parser = None

    @cached_method
    def one_basis(self):
        return (_coerce_index([1], False, self._ascode), self._varlist[0])

    @property
    def is_quantum(self):
        return self._quantum

    def set_coproduct_indices(self, indices):
        self._splitter = indices

    def _element_constructor_(self, x, vname=None):
        if isinstance(x, str):
            return self.parser().parse(x)
        if (x, vname) in self._index_wrapper:
            elem = self._from_dict({self._index_wrapper((x, vname)): self.base_ring().one()})
        elif isinstance(x, list) or isinstance(x, tuple) or isinstance(x, Permutation) or isinstance(x, Composition) or isinstance(x, pl.Permutation):
            elem = self._from_dict(
                {
                    self._index_wrapper(
                        (
                            _coerce_index(x, self._ascode, self._ascode),
                            self._varlist[0] if vname is None else vname,
                        ),
                    ): self.base_ring().one(),
                },
            )
        elif isinstance(x, FastDoubleSchubertPolynomial):
            if x.base_varname == self._base_varname and (self._quantum == x.parent()._quantum) and (not self._quantum or x.q_varname == self._q_varname):
                elem = self._from_dict(
                    {(_coerce_index(k[0], x.parent()._ascode, self._ascode), k[1]): v for k, v in x.monomial_coefficients().items()},
                )
            else:
                return self(x.expand())
        elif isinstance(x, bork.FastSchubertPolynomial):
            if x.base_varname == self._base_varname and x.q_varname == self._q_varname and (self._quantum == x.parent()._quantum) and (not self._quantum or x.q_varname == self._q_varname):
                elem_dict = {}

                for k, v in x.monomial_coefficients().items():
                    if self._quantum:
                        res = qyz.schubmult_q_double_fast(
                            {pl.Permutation([]): self._coeff_polynomial_ring(v)},
                            pl.Permutation(list(_coerce_index(k, x.parent()._ascode, False))),
                            self._coeff_polynomial_rings[self._varlist[0]].gens(),
                            [
                                0
                                for i in range(
                                    len(self._coeff_polynomial_rings[self._varlist[0]].gens()),
                                )
                            ],
                            self._q_ring.gens(),
                        )
                    else:
                        res = yz.schubmult_double(
                            {pl.Permutation([]): self._coeff_polynomial_ring(v)},
                            pl.Permutation(list(_coerce_index(k, x.parent()._ascode, False))),
                            self._coeff_polynomial_rings[self._varlist[0]].gens(),
                            [
                                0
                                for i in range(
                                    len(self._coeff_polynomial_rings[self._varlist[0]].gens()),
                                )
                            ],
                        )
                    for k0, c0 in res.items():
                        elem_dict[(_coerce_index(k0, False, self._ascode), self._varlist[0])] = elem_dict.get(
                            (_coerce_index(k0, False, self._ascode), self._varlist[0]),
                            self._coeff_polynomial_ring.zero(),
                        ) + self._coeff_polynomial_ring(c0)
                elem = self._from_dict(elem_dict)
            else:
                elem = self(x.expand())
        elif isinstance(x, MPolynomial):
            from sage.interfaces.sympy import sympy_init

            sympy_init()
            sympy_floff = sympify(str(x))
            val = syme.sympify(sympy_floff)
            if self._quantum:
                result = qyz.mult_poly_q_double(
                    {pl.Permutation([]): 1},
                    val,
                    [syme.Symbol(str(g)) for g in self._base_polynomial_ring.gens()],
                    [syme.Symbol(str(g)) for g in self._coeff_polynomial_rings[self._varlist[0]].gens()],
                    [syme.Symbol(str(g)) for g in self._q_ring.gens()],
                )
            else:
                result = yz.mult_poly_double(
                    {pl.Permutation([]): 1},
                    val,
                    [syme.Symbol(str(g)) for g in self._base_polynomial_ring.gens()],
                    [syme.Symbol(str(g)) for g in self._coeff_polynomial_rings[self._varlist[0]].gens()],
                )
            elem = self._from_dict(
                {
                    (
                        _coerce_index(k, False, self._ascode),
                        self._varlist[0],
                    ): self._coeff_polynomial_ring(str(v))
                    for k, v in result.items()
                },
            )
        else:
            raise TypeError

        return elem

    def some_elements(self):
        return [
            self.one(),
            self(_coerce_index([1], False, self._ascode)) + 2 * self(_coerce_index([2, 1], False, self._ascode)),
            self(_coerce_index([4, 2, 1, 3], False, self._ascode)) - self(_coerce_index([3, 2, 1], False, self._ascode)),
        ]

    def product_on_basis(self, left, right):
        le = _coerce_index(left[0], self._ascode, False)
        ri = _coerce_index(right[0], self._ascode, False)
        var_y = [syme.sympify(str(g)) for g in self._coeff_polynomial_rings[left[1]].gens()]
        var_z = [syme.sympify(str(g)) for g in self._coeff_polynomial_rings[right[1]].gens()]
        if self._quantum:
            result = qyz.schubmult_q_double_fast(
                {pl.Permutation(le): 1},
                pl.Permutation(ri),
                var_y,
                var_z,
                [syme.sympify(str(g)) for g in self._q_ring.gens()],
            )
        else:
            result = yz.schubmult_double(
                {pl.Permutation(le): 1},
                pl.Permutation(ri),
                var_y,
                var_z,
            )
        result = {k: v for k, v in result.items() if v != 0}
        return sum(
            [self._coeff_polynomial_ring(str(v)) * self(_coerce_index(k, False, self._ascode), left[1]) for k, v in result.items()],
        )

    def _coerce_map_from_(self, S):
        if isinstance(S, MPolynomialRing_base):
            return True
        if isinstance(S, bork.FastSchubertPolynomialRing_base):
            return True
        if isinstance(S, FastDoubleSchubertPolynomialRing_base):
            return True
        if isinstance(S, str):
            return True
        return super().has_coerce_map_from(S)

    def coproduct_on_basis(self, indm):
        if self._quantum:
            raise NotImplementedError("Quantum double Schubert polynomials have no coproduct")
        indices = self._splitter
        indices = sorted(indices)
        RR = self._coeff_polynomial_rings[indm[1]]
        TR = self._coeff_polynomial_rings[self._varlist[0]]

        coeff_dict = yz.schub_coprod_double(permtrim(indm[0]), indices, [syme.sympify(str(g)) for g in TR.gens()], [syme.sympify(str(g)) for g in RR.gens()])
        total_sum = 0
        for (fp, sp), val in coeff_dict.items():
            firstperm = Permutation(list(fp))
            secondperm = Permutation(list(sp))
            # val = syme.sympify(val).subs(subs_dict_coprod)
            total_sum += self._coeff_polynomial_ring(val) * self(_coerce_index(firstperm, False, self._ascode), indm[1]).tensor(self(_coerce_index(secondperm, False, self._ascode), self._varlist[0]))
        return total_sum

    def _repr_(self):
        return (
            f"Ring of {'quantum' if self._quantum else ''} double Schubert polynomials in {self._base_varname}{',' + self._q_varname if self._quantum else ''} with {len(self._base_polynomial_ring.gens())} variables with"  # noqa: E501
            f" coefficient variables {','.join(self._varlist)} over the ring {self._coeff_polynomial_ring.base_ring()} indexed by {'the Lehmer code' if self._ascode else 'permutations'}"
        )


FastDoubleSchubertPolynomial = FastDoubleSchubertPolynomial_class
FastDoubleSchubertPolynomialRing_base = FastDoubleSchubertPolynomialRing_xbasis
