"""Abstract base classes shared by every Schubert-family ring.

`BaseSchubertRing` holds a primary generating set (``genset``, the ``x`` variables)
and a coefficient generating set (``coeff_genset``, the ``y``/``z`` variables, or
``None`` for single Schubert polynomials), and defines the ring-level hooks that
concrete rings fill in: which multiplication kernel to use, how to expand a basis
element to a polynomial, how to print it, and how to change basis. `BaseSchubertElement`
is the corresponding dict-like element type (``{Permutation: coefficient}``).
"""

from sympy import pretty

from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import Add, CoercionFailed, S, expand, sympify, sympify_sympy, sympy_Mul
from schubmult.symbolic.common_polys import schubpoly_from_elems
from schubmult.utils._mul_utils import _mul_schub_dicts
from schubmult.utils.logging import get_logger
from schubmult.utils.perm_utils import add_perm_dict

from ..base_ring import BaseRing, BaseRingElement

logger = get_logger(__name__)


class BaseSchubertElement(BaseRingElement):
    """A linear combination of Schubert-family basis elements, stored as ``{Permutation: coeff}``."""

    def mult_poly(self, poly):
        """Multiply this element by an arbitrary polynomial ``poly`` in the ring's ``genset`` variables."""
        res_dict2 = {}
        for k, v in self.items():
            if self.ring.coeff_genset.label is None:
                dict2 = self.ring.mult_poly_single({k: v}, poly, self.ring.genset)
            else:
                dict2 = self.ring.mult_poly_double({k: v}, poly, self.ring.genset, self.ring.coeff_genset)
            res_dict2 = add_perm_dict(res_dict2, dict2)
        return self.ring.from_dict(res_dict2)

    def in_schubert_schur_basis(self, numvars):
        """Expand into the Schubert-tensor-Schur basis of the tensor square ring, splitting off the
        symmetric part in the last ``numvars`` variables.
        """
        res = (self.ring @ self.ring).zero
        for perm, v in self.items():
            res += v * self.ring.in_schubert_schur_basis(perm, numvars)
        return res

    def in_SEM_basis(self, elem_func=None):
        """Expand as a polynomial in elementary symmetric functions (the \"SEM\" presentation), using
        ``elem_func`` (default: the ring's symbolic ``symbol_elem_func``) as the elementary symmetric symbol.
        """
        if elem_func is None:
            elem_func = self.ring.symbol_elem_func
        result = S.Zero
        for k, v in self.items():
            result += sympify(v) * schubpoly_from_elems(k, self.ring.genset, self.ring.coeff_genset, elem_func=elem_func)
        return result

    def _sympystr(self, printer):
        return printer._print(pretty(self, use_unicode=False))

    def as_ordered_terms(self, *_, **__):
        """Terms ``coeff * basis_symbol`` sorted by permutation length then lexicographically (sympy printing hook)."""
        if len(self.keys()) == 0:
            return [sympify(S.Zero)]
        return [((self[k]) if k == self.ring.zero_monom else sympy_Mul(sympify_sympy(self[k]), self.ring.printing_term(k))) for k in sorted(self.keys(), key=(lambda kk: (kk.inv, tuple(kk)) if hasattr(kk, "inv") else kk))]

    def _eval_expand_basic(self, *args, **kwargs):  # noqa: ARG002
        return self.as_polynomial()

    def expand(self, deep=True, *args, **kwargs):  # noqa: ARG002
        """With ``deep=True`` (default) expand to an explicit polynomial in the variables; with
        ``deep=False`` only expand each coefficient, keeping the Schubert basis.
        """
        if not deep:
            return self.ring.from_dict({k: expand(v, **kwargs) for k, v in self.items()})
        return sympify(expand(self.as_polynomial()))

    def as_expr(self):
        """Sum of the ``as_terms()`` as a sympy ``Add``."""
        return Add(*self.as_terms())

    def as_polynomial(self):
        """Expand to an explicit polynomial: ``sum coeff * SchubertPoly(perm)``."""
        # print(f"{self=}")
        # try:
        return Add(*[v * self.ring.cached_schubpoly(k) for k, v in self.items()])
        # except SympifyError:
        #     return Add(*[sympify(v) * self.ring.cached_schubpoly(k) for k, v in self.items()])

    def as_classical(self):
        """Re-express in the classical (non-quantum) Schubert basis."""
        return self.ring.in_classical_basis(self)

    def as_quantum(self):
        """Re-express in the quantum Schubert basis."""
        return self.ring.in_quantum_basis(self)

    def almosteq(self, other):
        """Equality up to coefficient expansion (handles elements of different but compatible rings)."""
        if isinstance(other, BaseSchubertElement):
            elem1 = self
            elem2 = other
            if elem1.ring == elem2.ring:
                return (self - other).expand(deep=False).almosteq(S.Zero)
            return elem1.almosteq(elem1.ring.one * elem2)
        return (self - self.ring.from_expr(other)).expand(deep=False) == self.ring.zero

    def strip_zeros(self):
        """Drop basis elements whose coefficient is exactly zero."""
        return self.ring.from_dict({k: v for k, v in self.items() if v != S.Zero})


class BaseSchubertRing(BaseRing):
    """Abstract base ring for Schubert-family polynomials.

    Concrete subclasses supply the multiplication kernels (``double_mul``/``single_mul``,
    ``mult_poly_single``/``mult_poly_double``), the basis-element expansion
    (``cached_schubpoly``), printing (``printing_term``), coercion, and basis changes.
    Two rings compare equal iff they have the same type and generating sets.
    """

    def __eq__(self, other):
        return type(self) is type(other) and self.genset == other.genset and self.coeff_genset == other.coeff_genset

    def __init__(self, genset, coeff_genset, domain=None):
        """Args:
            genset: Primary generating set (the ``x`` variables).
            coeff_genset: Coefficient generating set (``y``/``z``), or a set with ``label=None`` for single rings.
            domain: Optional coefficient domain passed to `BaseRing`.
        """
        super().__init__(domain=domain)
        self._genset = genset
        self._coeff_genset = coeff_genset
        self.symbols = list(genset)
        self.zero_monom = Permutation([])

    def mul(self, elem, other):
        """Multiply two elements via `_mul_schub_dicts`, which dispatches to the appropriate
        `schubmult.mult` kernel based on both rings' generating sets; scalars go through `BaseRing.mul`.
        """
        if not isinstance(other, BaseSchubertElement):
            return super().mul(elem, other)
        prod = _mul_schub_dicts(elem, other, elem.ring, other.ring)
        # kernel output is already in the coefficient domain: skip the per-coefficient
        # free_symbols scan unless a subclass customizes from_dict
        if type(self).from_dict is BaseRing.from_dict:
            return self.from_dict_unchecked(prod)
        return self.from_dict(prod)


    def new(self, x):
        """Hook: build an element from ``x`` (permutation, code, or expression)."""

    def printing_term(self, k):
        """Hook: the sympy symbol displayed for basis element ``k``."""

    def coproduct_on_basis(self, k):
        """Hook: coproduct of basis element ``k`` in the tensor square ring."""

    def _coerce_mul(self, other):
        """Hook: coerce ``other`` for multiplication, or return ``None`` if not possible."""

    def is_elem_mul_type(self, elem):
        """Hook: whether ``elem`` should be multiplied via the elementary-symmetric fast path."""

    def elem_mul(self, ring_elem, elem):
        """Hook: elementary-symmetric fast-path multiplication."""

    def _coerce_add(self, other):
        """Hook: coerce ``other`` for addition, or return ``None`` if not possible."""

    @property
    def elem_sym(self):
        """Hook: the elementary symmetric polynomial function used by this ring."""

    @property
    def symbol_elem_func(self):
        """Hook: symbolic (unevaluated) elementary symmetric function for ``in_SEM_basis``."""

    def elem_sym_subs(self, kk):
        """Hook: substitution dict turning the symbolic elementary symmetric symbols back into polynomials."""

    def domain_new(self, element, orig_domain=None):  # noqa: ARG002
        """Coerce ``element`` into the coefficient domain, refusing anything containing a ``genset`` variable."""
        try:
            if isinstance(element, BaseRingElement):
                raise CoercionFailed("Not a domain element")
            element = sympify(element)
            if not any(x in self.genset for x in element.free_symbols):
                return element
            raise CoercionFailed("Ring element to coerce contains an element of the set of generators")
        except Exception:
            raise CoercionFailed(f"Could not coerce type {type(element)} to {self.__class__.__name__}")

    @property
    def genset(self):
        """Primary generating set (the ``x`` variables)."""
        return self._genset

    @property
    def coeff_genset(self):
        """Coefficient generating set (``y``/``z`` variables); ``label`` is ``None`` for single rings."""
        return self._coeff_genset

    def in_quantum_basis(self, elem):
        """Hook: re-express ``elem`` in the quantum Schubert basis."""

    def in_classical_basis(self, elem):
        """Hook: re-express ``elem`` in the classical Schubert basis."""

    def quantum_schubpoly(self, perm):
        """Hook: the quantum Schubert polynomial for ``perm``."""

    def cached_product(self, u, v, basis2):
        """Hook: cached structure constants of ``S_u * S_v`` with ``v`` in ``basis2``."""

    def cached_positive_product(self, u, v, basis2):
        """Hook: like ``cached_product`` but with manifestly positive coefficients."""

    def mul_expr(self, elem, x):
        """Hook: multiply ``elem`` by a symbolic expression ``x``."""

    @property
    def double_mul(self):
        """Hook: the double-Schubert multiplication kernel (e.g. ``schubmult_double``)."""

    @property
    def single_mul(self):
        """Hook: the single-Schubert multiplication kernel (e.g. ``schubmult_py``)."""

    @property
    def mult_poly_single(self):
        """Hook: the single-variant multiply-by-polynomial kernel."""

    @property
    def mult_poly_double(self):
        """Hook: the double-variant multiply-by-polynomial kernel."""

    @property
    def quantum_elem_func(self):
        """Hook: the quantum elementary symmetric function."""

    def cached_schubpoly(self, k):
        """Hook: the (cached) explicit polynomial for basis element ``k``."""
