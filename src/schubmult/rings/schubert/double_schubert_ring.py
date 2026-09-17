"""Double Schubert polynomial ring: the ``DSx`` interface.

`DoubleSchubertRing` represents ``Z[y][x]`` in the basis of double Schubert
polynomials ``S_w(x; y)``, dispatching products to `schubmult.mult.double`.
It is also the workhorse behind the single ring (`schubert_ring.SingleSchubertRing`
is a `DoubleSchubertRing` with an all-zero coefficient alphabet). Beyond ring
arithmetic, `DoubleSchubertElement` supports divided differences, isobaric
divided differences, variable substitution/evaluation, coproducts, and
expansion into elementary-symmetric (\"CEM\"/\"SEM\") bases.

Variants: `ElemDoubleSchubertRing` keeps coefficients as unevaluated factorial
elementary symmetric functions; `DoubleSchubertRingDown` uses the descent-side
(\"down\") kernels.
"""

from functools import cache, cached_property

import schubmult.mult.double as yz
import schubmult.mult.positivity as pos
import schubmult.mult.single as py
import schubmult.rings.printing as spolymod
import schubmult.utils.schub_lib as schub_lib
from schubmult.combinatorics.permutation import Permutation, uncode
from schubmult.symbolic import Add, DomainElement, Mul, Pow, S, Symbol, expand, expand_func, is_of_func_type, sympify, sympify_sympy
from schubmult.symbolic.common_polys import elem_sym_poly, schubpoly_classical_from_elems, schubpoly_from_elems, xreplace_genvars
from schubmult.symbolic.poly.variables import CustomGeneratingSet, GeneratingSet, GeneratingSet_base, MaskedGeneratingSet, NotEnoughGeneratorsError, genset_dict_from_expr, poly_genset
from schubmult.symbolic.symmetric_polynomials import CompleteSym_base, ElemSym_base, FactorialElemSym, coeffvars, degree, genvars, numvars, split_out_vars
from schubmult.utils.perm_utils import add_perm_dict

from ..tensor_ring import TensorRing
from .base_schubert_ring import BaseSchubertElement, BaseSchubertRing


def is_fact_elem_sym(obj):
    """Whether ``obj`` is an (unevaluated) factorial elementary symmetric function."""
    return is_of_func_type(obj, ElemSym_base)


def is_fact_complete_sym(obj):
    """Whether ``obj`` is an (unevaluated) factorial complete homogeneous symmetric function."""
    return is_of_func_type(obj, CompleteSym_base)


class DoubleSchubertElement(BaseSchubertElement):
    """An element of a `DoubleSchubertRing`: ``{Permutation: coefficient}`` in the
    double Schubert basis ``S_w(x; y)``, with sympy coefficients in ``y``.
    """

    def to_genset_dict(self, trim=False):
        """Expand to a polynomial and return ``{exponent_tuple: coeff}`` over the ``x`` variables;
        ``trim=True`` merges keys that differ only by trailing zeros.
        """
        gdict = genset_dict_from_expr(self.as_polynomial(), self.ring.genset)

        def _trim_tuple(tup):
            tup = [*tup]
            while len(tup) > 0 and tup[-1] == 0:
                tup.pop()
            return (*tup,)

        if trim:
            new_dict = {}
            for flop, val in gdict.items():
                real_tup = _trim_tuple(flop)
                new_dict[real_tup] = new_dict.get(real_tup, S.Zero) + val
        else:
            new_dict = gdict
        return new_dict

    def divdiff(self, i):
        """Divided difference ``partial_i``: ``S_w -> S_{w s_i}`` when ``i`` is a descent of ``w``, else 0."""
        return self.ring.from_dict({k.swap(i - 1, i): v for k, v in self.items() if i - 1 in k.descents()})

    def simpleref(self, i):
        """Action of the simple reflection ``s_i`` on the ``x`` variables: ``f + (x_{i+1} - x_i) partial_i f``."""
        return self + self.divdiff(i).mult_poly(self.ring.genset[i + 1] - self.ring.genset[i])

    def coeff_isobaric(self, i, beta):
        """Isobaric divided difference acting on the ``y`` (coefficient) alphabet, transported through
        the basis via the antipode-style inversion ``S_w -> (-1)^{l(w)} S_{w^{-1}}``.
        """
        coeff_ring = self.ring.coeff_ring
        result = self.ring.zero
        for k, v in self.items():
            coeff = coeff_ring.from_expr(v)
            coeff_term1 = coeff.simpleref(i).as_polynomial()
            coeff_term2 = coeff.isobaric(i, beta).as_polynomial() + beta * v

            basis_elemi = self.ring.from_dict({(~k): S.NegativeOne**k.inv})
            basis_term1i = basis_elemi.isobaric(i, beta) + beta * basis_elemi

            basis_term1 = self.ring.from_dict({(~k2): v2*(S.NegativeOne**k2.inv) for k2, v2 in basis_term1i.items()})
            basis_term2 = self.ring.from_dict({k: S.One})

            result += coeff_term1 * basis_term1 + coeff_term2 * basis_term2 - beta * v * self.ring(k)
        return result

    def isobaric(self, i, beta):
        """Beta-deformed isobaric divided difference ``pi_i = partial_i + beta (x_i partial_i - 1)``."""
        the_divdiff = self.divdiff(i)
        return the_divdiff + beta * (the_divdiff.mult_poly(self.ring.genset[i]) - self)

    def divdiff_perm(self, perm):
        """Apply ``partial_w`` for ``w = perm``, peeling simple reflections from the last descent."""
        if perm.inv == 0:
            return self
        desc = max(perm.descents())
        perm2 = perm.swap(desc, desc + 1)
        return self.divdiff(desc + 1).divdiff_perm(perm2)

    def isobaric_perm(self, perm, beta):
        """Apply the beta-isobaric ``pi_w`` for ``w = perm``."""
        if perm.inv == 0:
            return self
        desc = max(perm.descents())
        perm2 = perm.swap(desc, desc + 1)
        return self.isobaric_perm(perm2, beta).isobaric(desc + 1, beta)

    def isobaric_plus_beta(self, i, beta):
        """The variant ``partial_i + beta x_i partial_i`` (isobaric without the ``-beta`` identity term)."""
        the_divdiff = self.divdiff(i)
        return the_divdiff + beta * the_divdiff.mult_poly(self.ring.genset[i])

    def act(self, perm):
        """Permute the ``x`` variables by ``perm``, as a composition of ``simpleref``s."""
        perm = Permutation(perm)
        dset = perm.descents()
        if len(dset) == 0:
            return self
        i = next(iter(dset))
        return self.simpleref(i + 1).act(perm.swap(i, i + 1))

    def max_index(self):
        """The largest ``x`` index (1-indexed) any basis permutation actually depends on."""
        return max([max([0, *list(k.descents(zero_indexed=False))]) for k in self.keys()])

    def eval(self, x):
        """Substitute ``{generator: value}`` pairs one at a time (via ``pull_out_gen``); returns a
        scalar if the result collapses to the identity basis element.
        """
        ret = self
        for v, val in x.items():
            ret = ret.pull_out_gen(v)
            ret = ret.ring.from_dict({k: v2.subs(v, val) for k, v2 in ret.items()})
        if len(ret.keys()) == 1 and next(iter(ret.keys())) == Permutation([]):
            return ret[Permutation([])]
        return ret

    def subs(self, old, new):
        """Substitute ``old -> new`` where ``old`` is an ``x`` variable (moved to the last position and
        pulled out via ``pull_out_var``), a ``y`` variable (transported through the basis), or a plain
        coefficient symbol.
        """
        result = 0
        if self.ring.genset.index(old) != -1:
            result = 0
            index = self.ring.genset.index(old)
            mindex = self.max_index()
            if mindex < index:
                return self
            perm = Permutation([]).swap(index - 1, mindex)
            transf = self.act(perm)
            for k, v in transf.items():
                perm = k
                coeff_gens = self.ring.coeff_genset
                L = schub_lib.pull_out_var(mindex + 1, perm)
                for index_list, new_perm in L:
                    result += self.ring.from_dict({new_perm: v}).mult_poly(Mul(*[(new - coeff_gens[index2]) for index2 in index_list]))
            return result

        for k, v in self.items():
            if self.ring.coeff_genset.label is None:
                add_dict = {k: v.subs(old, new)}
            else:
                coeff_genset = self.ring.coeff_genset
                if coeff_genset.index(old) != -1:
                    genset_list = [coeff_genset[i] for i in range(len(coeff_genset))]
                    genset_list[coeff_genset.index(old)] = 0
                    custom_genset = CustomGeneratingSet(genset_list)
                    new_add_dict = {k2: sympify(v2).subs(old, new) for k2, v2 in yz.schubmult_double({(): v}, k, custom_genset, coeff_genset).items()}
                    add_dict = {}
                    for k3, v3 in new_add_dict.items():
                        to_add_dict = yz.schubmult_double({(): v3}, k3, coeff_genset, custom_genset)
                        add_dict = add_perm_dict(add_dict, to_add_dict)
                else:
                    add_dict = {k: sympify(v).subs(old, new)}
            for k5, v5 in add_dict.items():
                if any(self.ring.genset.index(s) != -1 for s in sympify(v5).free_symbols):
                    result += self.ring.from_dict({k5: 1}).mult_poly(v5)
                else:
                    result += self.ring.from_dict({k5: v5})
        return result

    @property
    def free_symbols(self):
        """Coefficient symbols plus the ``x``/``y`` variables the basis permutations actually depend on."""
        ret = set()
        for k, v in self.items():
            ret.update(v.free_symbols)
            perm = k
            if len(perm.descents()) > 0:
                ret.update([self.ring.genset[i] for i in range(1, max(perm.descents()) + 2)])
            if self.ring.coeff_genset:
                genset2 = self.ring.coeff_genset
                perm2 = ~perm
                if len(perm2.descents()) > 0:
                    ret.update([genset2[i] for i in range(1, max(perm2.descents()) + 2)])
        return ret

    def pull_out_gen(self, gen):
        """Factor out all dependence on one generator ``gen`` (an ``x`` or ``y`` variable), returning an
        element over a `MaskedGeneratingSet` ring with ``gen`` removed and explicit ``(gen - y_j)``
        (or factorial-elementary-symmetric) prefactors.
        """
        ind = self.ring.genset.index(gen)
        if ind == -1:
            ind = self.ring.coeff_genset.index(gen)
            if ind == -1:
                raise ValueError(f"{gen} passed but is not a generator")
            gens2 = MaskedGeneratingSet(self.ring.coeff_genset, [ind])
            gens2.set_label(f"({self.ring.coeff_genset.label}\\{gen})")
            new_basis = DoubleSchubertRing(self.ring.genset, gens2)
            ret = new_basis.zero
            for perm, val in self.items():
                L = schub_lib.pull_out_var(ind, ~perm)
                for index_list, new_perm in L:
                    ret += FactorialElemSym(len(index_list), len(index_list), [self.ring.genset[index2] for index2 in index_list], [gen]) * val * new_basis(~new_perm)
            return ret
        gens2 = MaskedGeneratingSet(self.ring.genset, [ind])
        gens2.set_label(f"({self.ring.genset.label}\\{gen})")
        new_basis = DoubleSchubertRing(gens2, self.ring.coeff_genset)
        ret = new_basis.zero
        for perm, val in self.items():
            L = schub_lib.pull_out_var(ind, perm)
            for index_list, new_perm in L:
                toadd = S.One
                for index2 in index_list:
                    toadd *= gen - self.ring.coeff_genset[index2]
                ret += toadd * val * new_basis(new_perm)
        return ret

    def in_CEM_basis(self):
        """Expand in the complete-elementary-monomial (CEM) basis using the ring's symbolic elementary function."""
        result = S.Zero
        for k, v in self.items():
            result += sympify(v) * schubpoly_classical_from_elems(k, self.ring.genset, self.ring.coeff_genset, elem_func=self.ring.symbol_elem_func)
        return result

    def cem_rep(self, elem_func, mumu=None):
        """CEM expansion with a custom ``elem_func``; ``mumu`` selects a dominant permutation to expand
        against (defaults to the classical route).
        """
        result = S.Zero
        if mumu is not None:
            for k, v in self.items():
                result += sympify(v) * schubpoly_from_elems(k, self.ring.genset, self.ring.coeff_genset, elem_func=elem_func, mumu=mumu)
        else:
            for k, v in self.items():
                result += sympify(v) * schubpoly_classical_from_elems(k, self.ring.genset, self.ring.coeff_genset, elem_func=elem_func)
        return result

    def coproduct(self, *indices, alt_coeff_genset=None, on_coeff_gens=False, gname1=None, gname2=None):
        """Coproduct splitting the ``x`` variables (or ``y`` if ``on_coeff_gens``) at the given 1-indexed
        ``indices``: returns an element of the `TensorRing` of two `DoubleSchubertRing`s over the
        complementary `MaskedGeneratingSet`s, labeled ``gname1``/``gname2``.
        """
        result_dict = {}
        genset = self.ring.genset
        if on_coeff_gens:
            genset = self.ring.coeff_genset
        if gname1 is None:
            gname1 = f"{genset.label}_A"
        if gname2 is None:
            gname2 = f"{genset.label}_B"
        gens2 = MaskedGeneratingSet(genset, indices)
        gens1 = gens2.complement()
        gens1.set_label(gname1)
        gens2.set_label(gname2)
        for k, v in self.items():
            key = k
            from .schubert_ring import SingleSchubertRing

            if isinstance(self.ring, SingleSchubertRing) and not alt_coeff_genset:
                coprod_dict = py.schub_coprod_py(key, indices)
            else:
                if on_coeff_gens:
                    coprod_dict = yz.schub_coprod_double(~key, indices, self.ring.genset, alt_coeff_genset if alt_coeff_genset else self.ring.genset)
                else:
                    coprod_dict = yz.schub_coprod_double(key, indices, self.ring.coeff_genset, alt_coeff_genset if alt_coeff_genset else self.ring.coeff_genset)
            if on_coeff_gens:
                result_dict = add_perm_dict(result_dict, {(~k1, ~k2): v * v2 for (k1, k2), v2 in coprod_dict.items()})
            else:
                result_dict = add_perm_dict(result_dict, {k: v * v2 for k, v2 in coprod_dict.items()})
        if on_coeff_gens:
            basis = TensorRing(
                DoubleSchubertRing(self.ring.genset, gens1),
                DoubleSchubertRing(alt_coeff_genset if alt_coeff_genset else self.ring.genset, gens2),
            )
        else:
            basis = TensorRing(
                DoubleSchubertRing(gens1, self.ring.coeff_genset),
                DoubleSchubertRing(gens2, alt_coeff_genset if alt_coeff_genset else self.ring.coeff_genset),
            )
        return basis.from_dict(result_dict)

    @cached_property
    def max_gens(self):
        """Largest 0-indexed descent over all basis permutations."""
        return max([max(k.descents()) for k in self.keys()])

    def positive_elem_sym_rep(self):
        """Manifestly positive expansion in factorial elementary symmetric functions (forward ``pull_out_var``)."""
        res = S.Zero
        for k, val in self.items():
            res += val * self.ring.positive_elem_sym_rep(k)
        return res

    def positive_elem_sym_rep_backward(self):
        """Like ``positive_elem_sym_rep`` but peeling from the last descent backward."""
        res = S.Zero
        for k, val in self.items():
            res += val * self.ring.positive_elem_sym_rep_backward(k)
        return res

    def antipode(self):
        """The antipode: swap the two alphabets and invert each basis permutation (see `DoubleSchubertRing.antipode`)."""
        return self.ring.antipode(self)


class DoubleSchubertRing(BaseSchubertRing):
    """The ring of double Schubert polynomials ``S_w(x; y)`` over ``genset`` (``x``) and
    ``coeff_genset`` (``y``). Call the ring with a permutation, Lehmer code, or polynomial
    expression to construct an element; the module-level ``DSx`` is the standard instance.
    """

    def __hash__(self):
        return hash((self.genset, self.coeff_genset, "DBS"))

    def __init__(self, genset, coeff_genset, domain=None):
        super().__init__(genset, coeff_genset, domain)
        self.dtype = type("DoubleSchubertElement", (DoubleSchubertElement,), {"ring": self})

    def __str__(self):
        return f"Double Schubert polynomial ring in {self.genset.label} and {self.coeff_genset.label}"

    @cached_property
    def coeff_ring(self):
        """The single Schubert ring over the coefficient alphabet ``y`` (used by ``coeff_isobaric``)."""
        from .schubert_ring import SingleSchubertRing
        return SingleSchubertRing(self.coeff_genset)

    @cached_property
    def antipode_ring(self):
        """The same ring with the two alphabets swapped."""
        return DoubleSchubertRing(self.coeff_genset, self.genset, domain=self.domain)

    def antipode(self, elem):
        """Map ``sum c_w S_w(x; y)`` to ``sum c_w S_{w^{-1}}(y; x)`` in the swapped-alphabet ring."""
        aring = self.antipode_ring
        result = aring.zero
        for k, v in elem.items():
            result += aring(~k) * v
        return result.expand(deep=False)

    def rmul(self, elem, other):
        """Right-multiply by a scalar (coefficient-domain element) or, failing that, by an expression."""
        try:
            other = self.domain_new(other)
            return self.from_dict({k: v * other for k, v in elem.items()})
        except Exception:
            return self.mul_expr(elem, other)

    def positive_elem_sym_rep(self, perm, index=1):
        """Manifestly positive expansion of ``S_perm`` in factorial elementary symmetric functions, peeling
        the first variable of ``~perm`` at each step (``pull_out_var(1, ...)``).
        """
        if perm.inv == 0:
            return S.One
        ret = S.Zero
        L = schub_lib.pull_out_var(1, ~perm)
        for index_list, new_perm in L:
            ret += self.elem_sym(len(index_list), len(index_list), [self.genset[index2] for index2 in index_list], [self.coeff_genset[index]]) * self.positive_elem_sym_rep(~new_perm, index + 1)
        return ret

    def positive_elem_sym_rep_backward(self, perm):
        """Like ``positive_elem_sym_rep`` but peeling from the last descent of ``~perm`` backward."""
        if perm.inv == 0:
            return S.One
        ret = S.Zero
        index = max((~perm).descents()) + 1
        L = schub_lib.pull_out_var(index, ~perm)
        for index_list, new_perm in L:
            ret += self.positive_elem_sym_rep_backward(~new_perm) * self.elem_sym(len(index_list), len(index_list), [self.genset[index2] for index2 in index_list], [self.coeff_genset[index]])
        return ret

    def printing_term(self, k, prefix=""):
        """The ``DSchubPoly`` display symbol for basis element ``k``."""
        return spolymod.DSchubPoly(k, self.genset.label, self.coeff_genset.label, prefix=prefix)

    def _coerce_mul(self, other):
        """Accept double/elem-double Schubert elements as-is; convert quantum double elements to classical."""
        from . import quantum_schubert_ring as qsr

        if isinstance(other, BaseSchubertElement):
            if isinstance(other.ring, qsr.QuantumDoubleSchubertRing):
                return other.as_classical()
            if isinstance(other.ring, ElemDoubleSchubertRing):
                return other
            if isinstance(other.ring, DoubleSchubertRing):
                return other
        return None

    def _coerce_add(self, other):  # noqa: ARG002
        return None

    @property
    def elem_sym(self):
        """`FactorialElemSym`."""
        return FactorialElemSym

    def is_elem_mul_type(self, other):
        """Whether ``other`` is a factorial elementary symmetric function (eligible for ``elem_mul``)."""
        return is_fact_elem_sym(other)

    def elem_mul(self, ring_elem, elem):
        """Multiply by a factorial elementary symmetric function in ``x`` variables via the positional
        Pieri rule (``elem_sym_positional_perms``), expanding the leftover factor with ``expand_func``.
        """
        elem = sympify(elem)
        indexes = [self.genset.index(a) for a in genvars(elem)]
        ret = self.zero
        elem_sympy = sympify_sympy(elem)
        for k, v in ring_elem.items():
            perm_list = schub_lib.elem_sym_positional_perms(k, degree(elem), *indexes)
            for perm, df, sign in perm_list:
                remaining_vars = [self.coeff_genset[perm[i - 1]] for i in indexes if perm[i - 1] == k[i - 1]]
                coeff = elem_sympy.func(degree(elem) - df, numvars(elem) - df, remaining_vars, coeffvars(elem))
                ret += (v * sign * expand_func(coeff)) * self(perm)
        return ret

    @property
    def symbol_elem_func(self):
        """`FactorialElemSym` (kept unevaluated for symbolic expansions)."""
        return FactorialElemSym

    def schubert_schur_elem_func(self, numvars):
        """Elementary-symmetric substitute for the Schubert-tensor-Schur expansion: ``e_p(x_1..x_k)`` maps
        to a Schubert basis element on the left factor when ``k >= numvars`` and on the right otherwise.
        """
        ring = self @ self

        def elem_func(p, k, *args):  # noqa: ARG001
            if p < 0:
                return ring.zero
            if p > k:
                return ring.zero
            if k >= 0 and p == 0:
                return ring.one
            if k >= numvars:
                return ring((uncode([0] * (k - p) + [1] * p), Permutation([])))
            return ring((Permutation([]), uncode([0] * (k - p) + [1] * p)))

        return elem_func

    def in_schubert_schur_basis(self, perm, numvars):
        """Expand ``S_perm`` in the Schubert-tensor-Schur basis, treating the last ``numvars`` variables
        as the symmetric (Schur) part.
        """
        elem_func = self.schubert_schur_elem_func(numvars)
        if perm.inv == 0:
            return elem_func(0, 0)
        extra = len(perm) - numvars
        dom = uncode([numvars] * extra + list(range(numvars - 1, 0, -1)))
        return schubpoly_from_elems(perm, self.genset, self.coeff_genset, elem_func=elem_func, mumu=dom)

    def in_descending_schur_basis(self, perm, numvars):
        """Iterate ``in_schubert_schur_basis`` down through ``numvars, numvars-1, ..., 1``, producing a
        nested tensor of Schur-like factors.
        """
        from .schubert_ring import Sx

        if numvars == 1:
            if perm == uncode([1]) or perm.inv == 0:
                return Sx(perm)
            return self.zero
        mid_res = self.in_schubert_schur_basis(perm, numvars)
        new_ring = TensorRing(mid_res.ring, self)
        result = new_ring.zero
        for k, v in mid_res.items():
            second_part = self.in_descending_schur_basis(k[1], numvars - 1)
            if numvars == 2:
                for k1, v2 in second_part.items():
                    result += new_ring.from_dict({(k[0], k1): v * v2})
            else:
                for (k1, k2), v2 in second_part.items():
                    result += new_ring.from_dict({((k[0], k1), k2): v * v2})
        return result

    def elem_sym_subs(self, kk):
        """Substitution dict ``{e_p_k: elem_sym_poly(p, k, x)}`` for all ``1 <= p <= k <= kk``."""
        elems = []
        for k in range(1, kk + 1):
            for p in range(1, k + 1):
                elems += [(Symbol(f"e_{p}_{k}"), elem_sym_poly(p, k, self.genset[1:], poly_genset(0)))]
        return dict(elems)

    @staticmethod
    def flip(elem):
        """Re-express a factorial elementary symmetric function with its two alphabets swapped, via the
        corresponding Grassmannian Schubert polynomial's CEM expansion.
        """
        R = DoubleSchubertRing(CustomGeneratingSet([0, *coeffvars(elem)]), CustomGeneratingSet([0, *genvars(elem)]))
        p = degree(elem)
        K = numvars(elem) + 1 - p
        poly = R(uncode([*list((K - 1) * [0]), p]))
        return poly.in_CEM_basis()

    def in_quantum_basis(self, elem):
        """Expand each basis element via ``quantum_schubpoly`` (a quantum double Schubert element)."""
        result = S.Zero
        for k, v in elem.items():
            result += v * self.quantum_schubpoly(k)
        return result

    def in_classical_basis(self, elem):
        """Identity (this ring is already classical)."""
        return elem

    @cache
    def quantum_schubpoly(self, perm):
        """The classical ``S_perm`` expressed in the quantum double Schubert basis (via ``quantum_elem_func``)."""
        return schubpoly_classical_from_elems(perm, self.genset, self.coeff_genset, self.quantum_elem_func)

    @cache
    def cached_product(self, u, v, basis2):
        """Structure constants of ``S_u(x; y) * S_v(x; z)`` (``z`` = ``basis2.coeff_genset``), via ``schubmult_double``."""
        return yz.schubmult_double({u: S.One}, v, self.coeff_genset, basis2.coeff_genset)

    @cache
    def cached_positive_product(self, u, v, basis2):
        """Like ``cached_product`` but with manifestly positive coefficients (generic alphabets, then substituted)."""
        return {k: xreplace_genvars(x, self.coeff_genset, basis2.coeff_genset) for k, x in pos.schubmult_generic_partial_posify(u, v).items()}

    @property
    def double_mul(self):
        """`schubmult.mult.double.schubmult_double`."""
        return yz.schubmult_double

    @property
    def single_mul(self):
        """`schubmult.mult.single.schubmult_py`."""
        return py.schubmult_py

    @property
    def mult_poly_single(self):
        """`schubmult.mult.single.mult_poly_py`."""
        return py.mult_poly_py

    @property
    def mult_poly_double(self):
        """`schubmult.mult.double.mult_poly_double`."""
        return yz.mult_poly_double

    @property
    def quantum_elem_func(self):
        """Elementary symmetric function valued in the quantum double Schubert ring, computed by a
        divide-and-conquer recursion on the variable set (used by ``quantum_schubpoly``).
        """
        from . import quantum_schubert_ring as qsr

        basis = qsr.QuantumDoubleSchubertRing(self.genset, self.coeff_genset)

        def elem_func(p, k, varl1, varl2, xstart=0, ystart=0):
            if p > k:
                return basis(0)
            if p == 0:
                return basis([])
            if p == 1:
                res = basis(varl1[xstart] - varl2[ystart])
                for i in range(1, k):
                    res += basis(varl1[xstart + i] - varl2[ystart + i])
                return res
            if p == k:
                res = basis((varl1[xstart] - varl2[ystart]) * (varl1[xstart + 1] - varl2[ystart]))
                for i in range(2, k):
                    res *= basis(varl1[i + xstart] - varl2[ystart])
                return res
            mid = k // 2
            xsm = xstart + mid
            ysm = ystart + mid
            kmm = k - mid
            res = elem_func(p, mid, varl1, varl2, xstart, ystart) + elem_func(
                p,
                kmm,
                varl1,
                varl2,
                xsm,
                ysm,
            )
            for p2 in range(max(1, p - kmm), min(p, mid + 1)):
                res += elem_func(p2, mid, varl1, varl2, xstart, ystart) * elem_func(
                    p - p2,
                    kmm,
                    varl1,
                    varl2,
                    xsm,
                    ysm - p2,
                )
            return res

        return elem_func

    def monomial_schub(self, monom):
        """The monomial ``x^monom`` expressed in the Schubert basis (trailing zeros in ``monom`` ignored)."""
        monom = [*monom]
        while len(monom) > 0 and monom[-1] == 0:
            monom.pop()
        return self._monomial_schub_cache(tuple(monom))

    @cache
    def _monomial_schub_cache(self, monom):
        """``monomial_schub`` worker: the dominant Schubert polynomial for the sorted exponent, permuted back."""
        srt_perm = Permutation.sorting_perm([-i for i in monom])
        schub_perm = uncode(sorted(monom, reverse=True))
        return self.from_dict({(schub_perm, 0): S.One}).act(srt_perm)

    @cache
    def cached_schubpoly(self, k):
        """The explicit polynomial ``S_k(x; y)`` (cached)."""
        return schubpoly_classical_from_elems(k, self.genset, self.coeff_genset, elem_func=elem_sym_poly)

    def complete_mul(self, elem, x):
        """Multiply by a factorial complete homogeneous symmetric function in ``x`` variables via
        ``complete_sym_positional_perms`` (the dual Pieri rule).
        """
        x = sympify(x)
        x_sympy = sympify_sympy(x)
        indexes = {self.genset.index(a) for a in genvars(x)}
        ret = self.zero
        for k, v in elem.items():
            perm_list = schub_lib.complete_sym_positional_perms(k, degree(x), *indexes)
            for perm, df, sign in perm_list:
                remaining_vars = [self.coeff_genset[perm[i - 1]] for i in {*indexes, *[j + 1 for j in range(len(perm)) if perm[j] != k[j]]}]
                coeff = x_sympy.func(degree(x) - df, numvars(x) + df, remaining_vars, coeffvars(x))
                ret += (sign * v * expand_func(coeff)) * self(perm)
        return ret

    def handle_sympoly(self, other):
        """How a symmetric-function coefficient is stored: evaluated to a polynomial here."""
        return expand_func(other)

    def single_variable(self, elem, varnum):
        """Multiply by the single variable ``x_varnum`` (equivariant Monk rule)."""
        ret = self.zero
        for u, v in elem.items():
            ret += v * self.coeff_genset[u[varnum - 1]] * self(u)
            new_perms = schub_lib.elem_sym_positional_perms(u, 1, varnum)
            for perm, udiff, sign in new_perms:
                if udiff == 1:
                    ret += (sign * v) * self(perm)
        return ret

    def from_expr(self, expr):
        """Convert a polynomial expression in ``x``/``y`` into the Schubert basis."""
        return super().from_expr(expr)

    def mul_expr(self, elem, x):
        """Multiply ``elem`` by an arbitrary expression ``x``: single variables use the Monk rule,
        (factorial) elementary/complete symmetric functions use their Pieri rules (splitting out
        variables from the wrong alphabet as needed), and ``Add``/``Mul``/``Pow`` recurse; anything
        else is treated as a coefficient.
        """
        if isinstance(x, DomainElement):
            raise TypeError(f"Cannot multiply {type(elem)} with {type(x)}")
        x = sympify(x)
        ind = self.genset.index(x)
        if ind != -1:
            return self.single_variable(elem, ind)
        if is_fact_elem_sym(x):
            if all(self.genset.index(a) != -1 for a in genvars(x)) and not any(self.genset.index(a) != -1 for a in coeffvars(x)):
                return self.elem_mul(elem, x)

            gens_to_remove = [a for a in genvars(x) if a not in self.genset]
            if any(self.genset.index(a) != -1 for a in genvars(x)) and len(gens_to_remove):
                return self.mul_expr(elem, x.split_out_vars(gens_to_remove))

            coeffs_to_remove = [a for a in coeffvars(x) if a in self.genset]

            if any(a in self.genset for a in coeffvars(x)) and len(coeffs_to_remove):
                return self.mul_expr(elem.split_out_vars(x.to_complete_sym(), coeffs_to_remove))
            return self.from_dict({k: (self.handle_sympoly(x)) * v for k, v in elem.items()})
        if is_fact_complete_sym(x):
            x = sympify(x)
            if all(a in self.genset for a in genvars(x)) and not any(a in self.genset for a in coeffvars(x)):
                return self.complete_mul(elem, x)
            gens_to_remove = [a for a in genvars(x) if a not in self.genset]

            if any(a in self.genset for a in genvars(x)) and len(gens_to_remove):
                return self.mul_expr(elem, x.split_out_vars(gens_to_remove))

            coeffs_to_remove = [a for a in x.coeff_vars if a in self.genset]

            if len(coeffs_to_remove):
                return self.mul_expr(elem, split_out_vars(x.to_elem_sym(), coeffs_to_remove))

            return self.from_dict({k: (self.handle_sympoly(x)) * v for k, v in elem.items()})
        if isinstance(x, Add):
            return self.sum([self.mul_expr(elem, arg) for arg in x.args])
        if isinstance(x, Mul):
            res = elem
            for arg in x.args:
                res = self.mul_expr(res, arg)
            return res
        if isinstance(x, Pow):
            res = elem
            base = x.args[0]
            exponent = x.args[1]
            if x.args[1] >= 0:
                #raise ValueError(f"Cannot multiply {elem} with {x}")
                for _ in range(int(exponent)):
                    res = self.mul_expr(res, base)
                return res
        #print(f"Fell through {x=}")
        return self.from_dict({k: v * self.domain_new(x) for k, v in elem.items()})

    def new(self, x):
        """Build an element from a permutation/Lehmer list, an existing element of this ring, or an expression."""
        genset = self.genset
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        if isinstance(x, list) or isinstance(x, tuple):
            p_x = Permutation(x)
            if max([0, *list(p_x.descents())]) > len(self.genset):
                raise NotEnoughGeneratorsError(f"Not enough generators {p_x=} {len(genset)=}")
            elem = self.from_dict({p_x: self.domain.one})
        elif isinstance(x, Permutation):
            if max([0, *list(x.descents())]) > len(self.genset):
                raise NotEnoughGeneratorsError(f"Not enough generators {x=} {len(genset)=}")
            elem = self.from_dict({x: self.domain.one})
        elif isinstance(x, DoubleSchubertElement):
            if x.ring.genset == genset:
                return x
            raise ValueError("Different generating set")
        else:
            elem = self.from_expr(x)
        return elem


class DoubleSchubertRingDown(DoubleSchubertRing):
    """`DoubleSchubertRing` using the descent-side (\"down\") multiplication kernels
    (``schubmult_double_down``/``schubmult_py_down``); basis symbols print with an ``op`` prefix.
    """

    def __hash__(self):
        return hash((self.genset, self.coeff_genset, "fatcabasi"))

    @property
    def double_mul(self):
        """`schubmult.mult.double.schubmult_double_down`."""
        return yz.schubmult_double_down

    @property
    def single_mul(self):
        """`schubmult.mult.single.schubmult_py_down`."""
        return py.schubmult_py_down

    @cache
    def cached_product(self, u, v, basis2):
        """Down-kernel structure constants over generic alphabets, substituted back to the ring's alphabets."""
        return {k: xreplace_genvars(x, self.coeff_genset, basis2.coeff_genset) for k, x in yz.schubmult_double_down({u: S.One}, v, yz._vars.var_g1, yz._vars.var_g2).items()}

    @cache
    def cached_positive_product(self, u, v, basis2):
        """Positive variant of ``cached_product`` for the down kernel."""
        return {k: xreplace_genvars(x, self.coeff_genset, basis2.coeff_genset) for k, x in pos.schubmult_double_down({u: S.One}, v, yz._vars.var_g1, yz._vars.var_g2).items()}

    def printing_term(self, k, prefix="op"):
        """The ``DSchubPoly`` display symbol, prefixed with ``op`` by default."""
        return spolymod.DSchubPoly(k, self.genset.label, self.coeff_genset.label, prefix=prefix)


class ElemDoubleSchubertRing(DoubleSchubertRing):
    """`DoubleSchubertRing` whose coefficients are kept as unevaluated `FactorialElemSym`
    functions instead of being expanded to polynomials; products use the ``*_from_elems`` kernels.
    """

    def __init__(self, genset, coeff_genset):
        super().__init__(genset, coeff_genset)
        self.dtype = type("DoubleSchubertElement", (DoubleSchubertElement,), {"ring": self})

    def __hash__(self):
        return hash((self.genset, self.coeff_genset, "EDBS"))

    @property
    def replacematch(self):
        """A ``(a, b) -> expression`` rewriter turning differences ``a - b`` into `FactorialElemSym(1, 1, ...)`
        forms, respecting which alphabet each symbol belongs to.
        """
        def bob(*args, **kwargs):  # noqa: ARG001
            a = kwargs["a"]
            b = kwargs["b"]
            ind1 = self.genset.index(a)
            ind2 = self.genset.index(b)
            if ind1 != -1:
                if ind2 == -1:
                    return FactorialElemSym(1, 1, [a], [b])
                return FactorialElemSym(1, 1, [a], [self.coeff_genset[1]]) - FactorialElemSym(1, 1, [b], [self.coeff_genset[1]])
            if ind2 != -1:
                return -FactorialElemSym(1, 1, [b], [a])
            if isinstance(a, Symbol) and isinstance(b, Symbol):
                return FactorialElemSym(1, 1, [a], [b])
            return a - b

        return bob

    @property
    def elem_func(self):
        """`FactorialElemSym`."""
        return FactorialElemSym

    def handle_sympoly(self, other):
        """Keep symmetric-function coefficients unevaluated."""
        return other

    def elem_mul(self, ring_elem, elem):
        """Positional Pieri rule for a factorial elementary symmetric function, keeping the leftover
        factor as an unevaluated coefficient.
        """
        indexes = [self.genset.index(a) for a in genvars(elem)]
        ret = self.zero
        elem_sympy = sympify_sympy(elem)
        for k, v in ring_elem.items():
            perm_list = schub_lib.elem_sym_positional_perms(k, degree(elem), *indexes)
            for perm, df, sign in perm_list:
                remaining_vars = [self.coeff_genset[perm[i - 1]] for i in indexes if perm[i - 1] == k[i - 1]]
                coeff = elem_sympy.func(degree(elem) - df, numvars(elem) - df, remaining_vars, coeffvars(elem))
                toadd = self.domain_new(v * sign * coeff) * self(perm)
                ret += toadd
        return ret

    def complete_mul(self, elem, x):
        """Dual Pieri rule for a factorial complete symmetric function, keeping the leftover factor unevaluated."""
        indexes = {self.genset.index(a) for a in genvars(x)}
        ret = self.zero
        x_sympy = sympify_sympy(x)
        for k, v in elem.items():
            perm_list = schub_lib.complete_sym_positional_perms(k, degree(x), *indexes)
            for perm, df, sign in perm_list:
                remaining_vars = [self.coeff_genset[perm[i - 1]] for i in {*indexes, *[j + 1 for j in range(len(perm)) if perm[j] != k[j]]}]
                coeff = x_sympy.func(degree(x) - df, numvars(x) + df, remaining_vars, coeffvars(x))
                ret += self.domain_new(sign * v * coeff) * self(perm)
        return ret

    @cache
    def cached_product(self, u, v, basis2):
        """Structure constants via ``schubmult_double_from_elems`` with `FactorialElemSym` coefficients."""
        return yz.schubmult_double_from_elems({u: self.domain.one}, v, self.coeff_genset, basis2.coeff_genset, elem_func=self.elem_func)

    @cache
    def cached_positive_product(self, u, v, basis2):
        """Structure constants via the positive ``schubmult_double_alt_from_elems`` route."""
        return {k: expand(v) for k, v in yz.schubmult_double_alt_from_elems({u: self.domain.one}, v, self.coeff_genset, basis2.coeff_genset, elem_func=self.elem_func).items()}

    def new(self, x):
        """Build an element from a permutation/Lehmer list, an element of this ring, or an expression."""
        genset = self.genset
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        if isinstance(x, list) or isinstance(x, tuple):
            p_x = Permutation(x)
            if max([0, *list(p_x.descents())]) > len(self.genset):
                raise NotEnoughGeneratorsError(f"Not enough generators {p_x=} {len(genset)=}")
            elem = self.from_dict({p_x: self.domain.one})
        elif isinstance(x, Permutation):
            if max([0, *list(x.descents())]) > len(self.genset):
                raise NotEnoughGeneratorsError(f"Not enough generators {x=} {len(genset)=}")
            elem = self.from_dict({x: self.domain.one})
        elif x.ring == self:
            return x
        else:
            elem = self.from_expr(x)
        return elem

    def _coerce_mul(self, other):
        """Accept double/elem-double Schubert elements as-is; convert quantum double elements to classical."""
        from . import quantum_schubert_ring as qsr

        if isinstance(other, BaseSchubertElement):
            if isinstance(other.ring, qsr.QuantumDoubleSchubertRing):
                return other.as_classical()
            if isinstance(other.ring, ElemDoubleSchubertRing):
                return other
            if isinstance(other.ring, DoubleSchubertRing):
                return other
        return None


def DSx(x, genset=GeneratingSet("y"), elem_sym=False, down=False):
    """Construct a double Schubert polynomial element in ``x`` with coefficient alphabet ``genset``.

    ``DSx([3, 1, 2])`` is ``S_{312}(x; y)``. Pass ``genset="z"`` (or a `GeneratingSet`) for a
    different coefficient alphabet; ``elem_sym=True`` uses `ElemDoubleSchubertRing`, ``down=True``
    uses `DoubleSchubertRingDown`.
    """
    if isinstance(genset, str):
        genset = GeneratingSet(genset)
    if down:
        return DoubleSchubertRingDown(GeneratingSet("x"), genset)(x)
    if elem_sym:
        return ElemDoubleSchubertRing(GeneratingSet("x"), genset)(x)
    return DoubleSchubertRing(GeneratingSet("x"), genset)(x)
