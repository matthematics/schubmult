"""Quantum double Schubert polynomial ring: the ``QDSx`` interface.

`QuantumDoubleSchubertRing` represents quantum double Schubert polynomials
``S^q_w(x; y)`` with quantum parameters ``q_1, q_2, ...`` (the module-level
``q_var``), dispatching products to `schubmult.mult.quantum_double`. Elements
can be converted to/from the classical basis via ``as_classical``/``as_quantum``.
"""

from functools import cache

import schubmult.mult.quantum as py
import schubmult.mult.quantum_double as yz
import schubmult.rings.schubert.schubert_ring as spr
import schubmult.utils.schub_lib as schub_lib
from schubmult.combinatorics.permutation import Permutation, uncode
from schubmult.symbolic import Add, Mul, Pow, S, Symbol, expand, expand_func, sympify
from schubmult.symbolic.common_polys import elem_sym_poly_q, schubpoly_from_elems, xreplace_genvars
from schubmult.symbolic.poly.variables import GeneratingSet, GeneratingSet_base, genset_dict_from_expr, poly_genset
from schubmult.utils._lazy import lazy_from

from ..printing import QDSchubPoly
from .base_schubert_ring import BaseSchubertElement, BaseSchubertRing

FactorialElemSym, QFactorialElemSym, coeffvars, degree, genvars, is_of_func_type, numvars = lazy_from(
    "schubmult.symbolic.symmetric_polynomials",
    "FactorialElemSym",
    "QFactorialElemSym",
    "coeffvars",
    "degree",
    "genvars",
    "is_of_func_type",
    "numvars",
)

q_var = GeneratingSet("q")


def is_fact_elem_sym(obj):
    """Whether ``obj`` is an (unevaluated) quantum factorial elementary symmetric function."""
    return is_of_func_type(obj, QFactorialElemSym)


class QuantumDoubleSchubertElement(BaseSchubertElement):
    """An element of a `QuantumDoubleSchubertRing`: ``{Permutation: coeff}`` in the basis ``S^q_w(x; y)``."""

    def subs(self, old, new):
        """Substitute by round-tripping through the classical basis."""
        return self.as_classical().subs(old, new).as_quantum()


class QuantumDoubleSchubertRing(BaseSchubertRing):
    """The ring of quantum double Schubert polynomials; ``QDSx`` is the standard constructor."""

    def __hash__(self):
        return hash((self.genset, self.coeff_genset, "QDBS"))

    def __init__(self, genset, coeff_genset):
        super().__init__(genset, coeff_genset)
        self.dtype = type("QuantumDoubleSchubertElement", (QuantumDoubleSchubertElement,), {"ring": self})

    def __str__(self):
        return f"Quantum Double Schubert polynomial ring in {self.genset.label} and {self.coeff_genset.label}"

    def printing_term(self, k):
        """The ``QDSchubPoly`` display symbol for basis element ``k``."""
        return QDSchubPoly(k, self.genset.label, self.coeff_genset.label)

    def _coerce_mul(self, other):
        """Accept quantum double elements over the same ``x``; lift quantum single elements by zeroing ``y``."""
        from .quantum_schubert_ring import QuantumSingleSchubertRing

        if isinstance(other, BaseSchubertElement):
            if other.ring == self:
                return other
            if isinstance(other.ring, QuantumDoubleSchubertRing) and other.ring.genset == self.genset:
                return other
            if type(other.ring) is QuantumSingleSchubertRing:
                if self.genset == other.ring.genset:
                    newbasis = QuantumDoubleSchubertRing(self.genset, poly_genset(0))
                    return newbasis.from_dict(other)
        return None

    def _coerce_add(self, other):
        """Accept only elements of an identical ring."""
        if isinstance(other, BaseSchubertElement):
            if type(other.ring) is type(self):
                if self.genset == other.ring.genset and self.coeff_genset == other.ring.coeff_genset:
                    return other
        return None

    @property
    def symbol_elem_func(self):
        """`QFactorialElemSym`."""
        return QFactorialElemSym

    def elem_sym_subs(self, kk):
        """Substitution dict ``{e_p_k: elem_sym_poly_q(p, k, x)}`` for all ``1 <= p <= k <= kk``."""
        elems = []
        for k in range(1, kk + 1):
            for p in range(1, k + 1):
                elems += [(Symbol(f"e_{p}_{k}"), elem_sym_poly_q(p, k, self.genset[1:], poly_genset(0)))]
        return dict(elems)

    @property
    def elem_sym(self):
        """`QFactorialElemSym`."""
        return QFactorialElemSym

    def is_elem_mul_type(self, other):
        """Whether ``other`` is a quantum factorial elementary symmetric function."""
        return is_fact_elem_sym(other)

    def elem_mul(self, ring_elem, elem):
        """Multiply by a quantum factorial elementary symmetric function via the quantum positional
        Pieri rule (``elem_sym_positional_perms_q``), each term carrying its ``q``-monomial.
        """
        elem = sympify(elem)
        indexes = [self.genset.index(a) for a in genvars(elem)]
        ret = self.zero
        for k, v in ring_elem.items():
            perm_list = schub_lib.elem_sym_positional_perms_q(k, degree(elem), *indexes, q_var=q_var)
            for perm, df, sign, mul_val in perm_list:
                remaining_vars = [self.coeff_genset[perm[i - 1]] for i in indexes if perm[i - 1] == k[i - 1]]
                coeff = FactorialElemSym(degree(elem) - df, numvars(elem) - df, remaining_vars, coeffvars(elem))
                ret += (v * sign * expand_func(coeff) * mul_val) * self(perm)
        return ret

    @cache
    def cached_product(self, u, v, basis2):
        """Structure constants via ``schubmult_q_double_fast``."""
        return yz.schubmult_q_double_fast({u: S.One}, v, self.coeff_genset, basis2.coeff_genset)

    def in_quantum_basis(self, elem):
        """Identity (this ring is already quantum)."""
        return elem

    def in_classical_basis(self, elem):
        """Expand each ``S^q_w`` as a classical double Schubert element via ``quantum_as_classical_schubpoly``."""
        result = None
        for k, v in elem.items():
            if not result:
                result = v * self.quantum_as_classical_schubpoly(k)
            else:
                result += v * self.quantum_as_classical_schubpoly(k)
        return result if result else self.zero

    @property
    def classical_elem_func(self):
        """Quantum elementary symmetric function valued in the classical `DoubleSchubertRing`, via the
        recursion ``E_p(k) = (x_k - y_{k-p+1}) E_{p-1}(k-1) + E_p(k-1) + q_{k-1} E_{p-2}(k-2)``.
        """
        basis = spr.DoubleSchubertRing(self.genset, self.coeff_genset)
        qv = yz._vars.q_var

        def elem_func(p, k, varl1, varl2):
            if p == 0 and k >= 0:
                return basis.one
            if p < 0 or p > k:
                return basis.zero
            return (varl1[k - 1] - varl2[k - p]) * elem_func(p - 1, k - 1, varl1, varl2) + elem_func(p, k - 1, varl1, varl2) + qv[k - 1] * elem_func(p - 2, k - 2, varl1, varl2)

        return elem_func

    @cache
    def quantum_as_classical_schubpoly(self, perm):
        """``S^q_perm`` expanded in the classical double Schubert basis (cached)."""
        return schubpoly_from_elems(perm, self.genset, self.coeff_genset, self.classical_elem_func)

    @cache
    def cached_schubpoly(self, k):
        """The explicit quantum double Schubert polynomial ``S^q_k(x; y)`` (cached)."""
        return schubpoly_from_elems(k, self.genset, self.coeff_genset, elem_func=elem_sym_poly_q)

    @cache
    def cached_positive_product(self, u, v, basis2):
        """Structure constants with (partially) positive coefficients via ``schubmult_q_generic_partial_posify``."""
        return {k: xreplace_genvars(x, self.coeff_genset, basis2.coeff_genset) for k, x in yz.schubmult_q_generic_partial_posify(u, v).items()}

    @property
    def double_mul(self):
        """`schubmult.mult.quantum_double.schubmult_q_double_fast`."""
        return yz.schubmult_q_double_fast

    @property
    def single_mul(self):
        """`schubmult.mult.quantum.schubmult_q_fast`."""
        return py.schubmult_q_fast

    @property
    def mult_poly_single(self):
        """`schubmult.mult.quantum.mult_poly_q`."""
        return py.mult_poly_q

    def positive_elem_sym_rep(self, perm, index=1):
        """Manifestly positive expansion of ``S^q_perm`` in quantum factorial elementary symmetric functions (forward)."""
        if perm.inv == 0:
            return S.One
        ret = S.Zero
        L = schub_lib.pull_out_var(1, ~perm)
        for index_list, new_perm in L:
            ret += self.elem_sym(len(index_list), len(index_list), [self.genset[index2] for index2 in index_list], [self.coeff_genset[index]]) * self.positive_elem_sym_rep(~new_perm, index + 1)
        return ret

    def positive_elem_sym_rep_backward(self, perm):
        """Like ``positive_elem_sym_rep`` but peeling from the last descent backward."""
        if perm.inv == 0:
            return S.One
        ret = S.Zero
        index = max((~perm).descents()) + 1
        L = schub_lib.pull_out_var(index, ~perm)
        for index_list, new_perm in L:
            ret += self.positive_elem_sym_rep_backward(~new_perm) * self.elem_sym(len(index_list), len(index_list), [self.genset[index2] for index2 in index_list], [self.coeff_genset[index]])
        return ret

    @property
    def mult_poly_double(self):
        """`schubmult.mult.quantum_double.mult_poly_q_double`."""
        return yz.mult_poly_q_double

    def from_expr(self, expr):
        """Convert a polynomial to the quantum basis by repeatedly peeling off the leading monomial's
        Schubert term; falls back to ``mul_expr`` on the identity if that fails.
        """
        ret = self.zero
        try:
            expr = sympify(expr)
            while expr != S.Zero:
                dct = genset_dict_from_expr(expr, self.genset)
                key = sorted(dct.keys(), reverse=True)[0]
                term = self.from_dict({uncode(key): dct[key]})
                ret += term
                expr -= term.as_polynomial()
                expr = expand(expr)
            return ret
        except Exception:
            return self.mul_expr(self.one, expr)

    def handle_sympoly(self, other):
        """Evaluate symmetric-function coefficients to polynomials."""
        return expand_func(other)

    def mul_expr(self, elem, x):
        """Multiply by an expression: single ``x`` variables via ``mult_poly_q_double``, quantum factorial
        elementary symmetric functions via ``elem_mul``, ``Add``/``Mul``/``Pow`` recursively, else as a coefficient.
        """
        x = sympify(x)
        _Add = Add
        _Mul = Mul
        _Pow = Pow
        ind = self.genset.index(x)
        if ind != -1:
            return self.from_dict(yz.mult_poly_q_double(elem, x, self.genset, self.coeff_genset))
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
        if isinstance(x, _Add):
            return self.sum([self.mul_expr(elem, arg) for arg in x.args])
        if isinstance(x, _Mul):
            res = elem
            for arg in x.args:
                res = self.mul_expr(res, arg)
            return res
        if isinstance(x, _Pow):
            res = elem
            for _ in range(int(x.args[1])):
                res = self.mul_expr(res, x.args[0])
            return res
        return self.from_dict({k: v * self.domain_new(x) for k, v in elem.items()})

    def new(self, x):
        """Build an element from a permutation/Lehmer list or a polynomial expression."""
        genset = self.genset
        if not genset:
            genset = self.genset
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        if isinstance(x, list) or isinstance(x, tuple):
            elem = self.from_dict({Permutation(x): self.domain.one})
        elif isinstance(x, Permutation):
            elem = self.from_dict({x: self.domain.one})
        else:
            elem = self.from_expr(x)
        return elem


def QDSx(x, genset=GeneratingSet("y")):
    """Construct a quantum double Schubert element in ``x`` with coefficient alphabet ``genset``
    (a `GeneratingSet` or a label string); e.g. ``QDSx([3, 1, 2])`` is ``S^q_{312}(x; y)``.
    """
    if isinstance(genset, str):
        genset = GeneratingSet(genset)
    return QuantumDoubleSchubertRing(GeneratingSet("x"), genset)(x)


QuantumDoubleSchubertPolynomial = QuantumDoubleSchubertElement
