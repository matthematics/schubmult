"""Shared machinery for every ring in `schubmult.rings`.

`BaseRingElement` is a ``dict`` mapping basis keys (permutations, RC graphs, tuples,
...) to coefficients, wired into sympy's printing and arithmetic protocols so that
elements can be added, multiplied, and displayed. `BaseRing` provides the
corresponding ring-level operations (``add``/``sub``/``mul``, coercion via
``domain_new``, construction via ``from_dict``/``from_expr``) and declares the hooks
concrete rings must implement (``new``, ``printing_term``, ``mul_expr``, ...).
A ring's element type is created dynamically as ``self.dtype`` with ``ring`` bound.
"""

from schubmult.symbolic import EXRAW, CoercionFailed, CompositeDomain, DefaultPrinting, DomainElement, Ring, S, expand, sympify, sympify_sympy, sympy_Add, sympy_Mul
from schubmult.utils.perm_utils import add_perm_dict


class BaseRingElement(DomainElement, DefaultPrinting, dict):
    """A ring element: ``{basis_key: coefficient}`` with sympy-compatible arithmetic and printing."""

    _op_priority = 1e200
    precedence = 40

    __sympy__ = True

    def __reduce__(self):
        return (self.__class__, self.items())

    @property
    def is_zero(self):
        """Whether every coefficient is exactly zero."""
        return all(v == S.Zero for v in self.values())

    def parent(self):
        """The ring this element belongs to (sympy domain protocol)."""
        return self.ring

    def has_free(self, *args):
        """Whether any of the given symbols appears in ``free_symbols``."""
        return any(s in args for s in self.free_symbols)

    def apply_to_keys(self, func):
        """Map each basis key through ``func`` (dropping keys where it returns ``None``), keeping coefficients."""
        new_elem = self.ring.zero
        for k, v in self.items():
            new_k = func(k)
            if new_k is None:
                continue
            new_elem += v * self.ring(new_k)
        return new_elem

    def eval(self, *args):
        pass

    def _sympystr(self, printer):
        if len(self.keys()) == 0:
            return printer._print(S.Zero)
        if printer.order in ("old", "none"):
            return printer._print_Add(sympy_Add(*self.as_ordered_terms()), order="lex")
        return printer._print_Add(sympy_Add(*self.as_ordered_terms()))

    def _pretty(self, printer):
        if len(self.keys()) == 0:
            return printer._print(S.Zero)
        if printer.order in ("old", "none"):
            return printer._print_Add(self, order="lex")
        return printer._print_Add(sympy_Add(*self.as_ordered_terms()))

    def _latex(self, printer):
        if len(self.keys()) == 0:
            return printer._print(S.Zero)
        if printer.order in ("old", "none"):
            return printer._print_Add(self, order="lex")
        return printer._print_Add(sympy_Add(*self.as_ordered_terms()))

    def as_terms(self):
        """Terms ``coeff * basis_symbol`` in dict order (sympy printing hook)."""
        if len(self.keys()) == 0:
            return [sympify_sympy(S.Zero)]
        return [((self[k]) if k == self.ring.zero_monom else sympy_Mul(sympify_sympy(self[k]), self.ring.printing_term(k))) for k in self.keys()]

    def as_ordered_terms(self, *_, **__):
        """Terms sorted by basis key (sympy printing hook)."""
        if len(self.keys()) == 0:
            return [sympify(S.Zero)]
        try:
            return [((self[k]) if k == self.ring.zero_monom else sympy_Mul(sympify_sympy(self[k]), self.ring.printing_term(k))) for k in sorted(self.keys())]
        except Exception as e:
            print([self[k] for k in sorted(self.keys())])
            raise e

    def __add__(self, other):
        if isinstance(other, BaseRingElement):
            return self.ring.add(self, other)
        try:
            other = self.ring.domain_new(other)
            other = self.ring.from_dict({self.ring.zero_monom: other})
            return self.ring.add(self, other)
        except CoercionFailed:
            pass
        try:
            new_other = self.ring(other)
            return self.__add__(new_other)
        except CoercionFailed:
            return other.__radd__(self)

    def __radd__(self, other):
        try:
            other = self.ring.domain_new(other)
            other = self.ring.from_dict({self.ring.zero_monom: other})
            return self.ring.add(other, self)
        except CoercionFailed:
            pass
        try:
            new_other = self.ring(other)
            return new_other.__add__(self)
        except CoercionFailed:
            return NotImplemented

    def __sub__(self, other):
        if isinstance(other, BaseRingElement):
            return self.ring.sub(self, other)
        try:
            other = self.ring.domain_new(other)
            other = self.ring.from_dict({self.ring.zero_monom: other})
            return self.ring.sub(self, other)
        except CoercionFailed:
            pass
        try:
            new_other = self.ring(other)
            return self.__sub__(new_other)
        except CoercionFailed:
            return other.__rsub__(self)

    def __rsub__(self, other):
        try:
            other = self.ring.domain_new(other)
            other = self.ring.from_dict({self.ring.zero_monom: other})
            return self.ring.sub(other, self)
        except CoercionFailed:
            return NotImplemented

    def __neg__(self):
        return self.ring.neg(self)

    def __pow__(self, val):
        try:
            val = int(val)
        except Exception:
            return NotImplemented
        if val == 0:
            return self.ring.one
        if val < 0:
            return NotImplemented
        if val == 1:
            return self
        return (self ** (val - 1)) * self

    def __mul__(self, other):
        return self.ring.mul(self, other)

    def __rmul__(self, other):
        try:
            return self.ring.rmul(self, other)
        except CoercionFailed:
            return NotImplemented

    def coproduct(self):
        """Coproduct into the tensor square ring, via ``ring.coproduct_on_basis``."""
        result = (self.ring@self.ring).zero
        for k, v in self.items():
            result += v * self.ring.coproduct_on_basis(k)
        return result

    def as_coefficients_dict(self):
        """``{basis_symbol: coeff}`` mapping display symbols to coefficients."""
        return {self.ring.printing_term(k, self.ring): sympify(v) for k, v in self.items()}

    def _eval_expand_basic(self, *args, **kwargs):  # noqa: ARG002
        return self.as_polynomial()

    def expand(self, deep=True, *args, **kwargs):  # noqa: ARG002
        """With ``deep=True`` expand to an explicit polynomial (``as_polynomial``); with ``deep=False`` only
        expand each coefficient, keeping the basis.
        """
        if not deep:
            return self.ring.from_dict({k: expand(v, **kwargs) for k, v in self.items()})
        return sympify(expand(self.as_polynomial()))

    def as_expr(self):
        """Sum of the ``as_terms()`` as a sympy ``Add``."""
        return sympy_Add(*self.as_terms())

    def as_polynomial(self):
        """Hook: expand this element to an explicit polynomial expression."""


    def __eq__(self, other):
        return type(self) is type(other) and self.ring == other.ring and dict.__eq__(self, other)

    def almosteq(self, other):
        """Equality up to coefficient expansion."""
        if isinstance(other, BaseRingElement):
            elem1 = self
            elem2 = other
            test_elem = elem1 - elem2
            if all(expand(v) == 0 for v in test_elem.values()):
                return True
            return False
        return False

    def __matmul__(self, other):
        """Tensor product ``self (x) other`` in the `TensorRing` of the two rings."""
        return (self.ring @ other.ring).ext_multiply(self, other)


class BaseRing(Ring, CompositeDomain):
    """Abstract base ring over a sympy coefficient domain (default ``EXRAW``); see the module docstring."""

    def __str__(self):
        return self.__class__.__name__

    def __matmul__(self, other):
        """The `TensorRing` ``self (x) other``."""
        from .tensor_ring import TensorRing

        return TensorRing(self, other)

    def __eq__(self, other):
        return type(self) is type(other)

    def to_sympy(self, elem):
        """Convert an element to a sympy expression (``as_expr``)."""
        return elem.as_expr()

    def __init__(self, domain=None):
        """Args:
            domain: Coefficient domain; defaults to sympy's ``EXRAW`` (arbitrary expressions).
        """
        if domain:
            self.domain = domain
        else:
            self.domain = EXRAW
        self.dom = self.domain
        self.zero_monom = ()

    def add(self, elem, other):
        """Coefficient-wise sum, dropping zeros."""
        res = self.from_dict(add_perm_dict(elem, other))
        return self.from_dict({k: v for k, v in res.items() if v != S.Zero})

    def sub(self, elem, other):
        """Coefficient-wise difference, dropping zeros."""
        res = self.from_dict(add_perm_dict(elem, {k: -v for k, v in other.items()}))
        return self.from_dict({k: v for k, v in res.items() if v != S.Zero})

    def neg(self, elem):
        """Negate every coefficient."""
        return self.from_dict({k: -v for k, v in elem.items()})

    def rmul(self, elem, other):
        """Right-multiply by a scalar (via ``domain_new``), falling back to ``mul_expr``."""
        try:
            other = self.domain_new(other)
            return self.from_dict({k: v * other for k, v in elem.items()})
        except CoercionFailed:
            return self.mul_expr(elem, other)

    def mul(self, elem, other):
        """Multiply by a scalar (via ``domain_new``), falling back to ``mul_expr``."""
        try:
            other = self.domain_new(other)
            return self.from_dict({k: other * v for k, v in elem.items()})
        except CoercionFailed:
            return self.mul_expr(elem, other)

    def to_domain(self):
        return self

    def from_sympy(self, expr):
        """Alias for ``from_expr`` (sympy domain protocol)."""
        return self.from_expr(expr)

    def new(self, x):
        """Hook: build an element from ``x``."""

    def printing_term(self, k):
        """Hook: the sympy symbol displayed for basis key ``k``."""

    def coproduct_on_basis(self, k):
        """Hook: coproduct of basis key ``k`` in the tensor square ring."""

    def _coerce_mul(self, other):  # noqa: ARG002
        """Hook: coerce ``other`` for multiplication, or ``None``."""
        return

    @property
    def one(self):
        """The multiplicative identity: coefficient 1 on ``zero_monom``."""
        return self.from_dict({self.zero_monom: S.One})

    def is_elem_mul_type(self, elem):
        """Hook: whether ``elem`` should use the ``elem_mul`` fast path."""

    def elem_mul(self, ring_elem, elem):
        """Hook: elementary-symmetric fast-path multiplication."""

    def _coerce_add(self, other):  # noqa: ARG002
        """Hook: coerce ``other`` for addition, or ``None``."""
        return

    def from_dict(self, element, orig_domain=None):
        """Build an element from ``{key: coeff}``, coercing each coefficient via ``domain_new`` and dropping zeros."""
        domain_new = self.domain_new
        poly = self.zero

        for monom, coeff in element.items():
            coeff = domain_new(coeff, orig_domain)
            if coeff != self.domain.zero:
                poly[monom] = coeff
        return poly

    def from_dict_unchecked(self, element):
        """from_dict for coefficients already known to lie in the domain (drops structural zeros only)."""
        poly = self.zero
        zero = self.domain.zero
        for monom, coeff in element.items():
            if coeff != zero:
                poly[monom] = coeff
        return poly

    @property
    def zero(self):
        """The empty element."""
        return self.dtype()

    def domain_new(self, element, orig_domain=None):  # noqa: ARG002
        """Coerce ``element`` into the coefficient domain (``sympify``), refusing ring/domain elements."""
        try:
            if isinstance(element, BaseRingElement) or isinstance(element, DomainElement):
                raise CoercionFailed("Not a domain element")
            return sympify(element)
        except Exception:
            raise CoercionFailed(f"Could not coerce type {type(element)} to {self.__class__.__name__}")

    def from_expr(self, x):
        """Build an element from an expression by multiplying the identity by it."""
        return self.mul_expr(self.one, x)

    def mul_expr(self, elem, x):
        """Hook: multiply ``elem`` by a symbolic expression ``x``."""


