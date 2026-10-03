r"""
Unexpanded coefficients: a Sage ring of raw SymEngine expressions

Every native Sage polynomial ring keeps its elements in expanded normal form. For Schubert structure
constants that is an exponential loss: the kernels produce coefficients as products of linear
factors (``(y_1 - y_3)*(y_2 - y_3)``, ``prod (1 + beta*y_i)**e`` denominators) whose expansions are
vastly larger. :class:`SymEngineRing` is a coefficient ring whose elements are SymEngine expression
trees kept exactly as the kernel built them: ``+`` and ``*`` are structural, zero is the literal
``0``, and nothing is ever expanded unless asked for. A
:class:`~sage.combinat.free_module.CombinatorialFreeModule` over it removes zero terms by truth
value, so raw kernel output can be used as module coefficients without any conversion.

Variables keep schubmult's names and 1-based indices (``y_1`` is the first variable), unlike the
polynomial base rings of the other parents in :mod:`schubmult.sage`, which are 0-indexed.

EXAMPLES::

    sage: from schubmult.sage.symengine_ring import SymEngineRing
    sage: E = SymEngineRing(QQ); E
    Ring of SymEngine expressions over Rational Field
    sage: y1, y2, y3 = (E.variable('y', i) for i in (1, 2, 3))
    sage: c = (y1 - y3) * (y2 - y3); c
    (y_1 - y_3)*(y_2 - y_3)
    sage: c.expand()
    y_2*y_1 - y_2*y_3 - y_3*y_1 + y_3**2
    sage: c == c.expand()                          # equality of two expressions is semantic ...
    True
    sage: bool(c - c.expand()), c - c.expand() == 0   # ... zero tests are structural
    (True, False)
    sage: M = CombinatorialFreeModule(E, Permutations())
    sage: v = c * M.monomial(Permutation([2, 1])) + (y1 - y2) * M.monomial(Permutation([1, 3, 2])); v
    (y_1-y_2)*B[[1, 3, 2]] + ((y_1-y_3)*(y_2-y_3))*B[[2, 1]]
    sage: v.map_coefficients(lambda a: a.expand())
    (y_1-y_2)*B[[1, 3, 2]] + (y_2*y_1-y_2*y_3-y_3*y_1+y_3**2)*B[[2, 1]]
    sage: c.to_polynomial(InfinitePolynomialRing(QQ, ['y']))
    y_2^2 - y_2*y_1 - y_2*y_0 + y_1*y_0

Expressions can be moved into a Sage polynomial ring (there the indices become 0-based), and Sage
polynomials and rationals coerce in::

    sage: R.<a0, a1> = QQ[]
    sage: E(a0*a1 + 1) * c
    (y_1 - y_3)*(y_2 - y_3)*(1 + a_2*a_1)
    sage: E(QQ(1)/2) * y1
    (1/2)*y_1
"""

from fractions import Fraction
from random import Random

from sage.categories.commutative_rings import CommutativeRings
from sage.misc.cachefunc import cached_method
from sage.rings.fraction_field import FractionField_generic
from sage.rings.fraction_field_element import FractionFieldElement
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.infinite_polynomial_element import InfinitePolynomial
from sage.rings.polynomial.infinite_polynomial_ring import InfinitePolynomialRing_sparse
from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.polynomial.multi_polynomial_ring_base import MPolynomialRing_base
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.polynomial.polynomial_ring import PolynomialRing_generic
from sage.rings.rational import Rational
from sage.rings.rational_field import QQ
from sage.structure.element import CommutativeRingElement
from sage.structure.parent import Parent
from sage.structure.richcmp import op_EQ, op_NE
from sage.structure.unique_representation import UniqueRepresentation

from ._convert import parse_sage_name, sage_polynomial_to_symengine, symengine_to_base_ring

_TRIAL_POINTS = 3
_TRIAL_RANGE = 10**9
_random = Random(0x5EED)  # noqa: S311  (reproducibility, not security)


def is_identically_zero(expr):
    """Whether the SymEngine expression ``expr`` is zero as a rational function.

    Evaluation at a few random integer points first (linear in the size of the tree): a nonzero
    value proves nonvanishing, which is the common case. Only an expression vanishing at every
    point is expanded (over a common denominator) for a proof.
    """
    import symengine

    if expr == 0:
        return True
    symbols = list(expr.free_symbols)
    for _ in range(_TRIAL_POINTS):
        value = expr.subs({s: symengine.Integer(_random.randrange(-_TRIAL_RANGE, _TRIAL_RANGE)) for s in symbols})
        if value.is_Integer or value.is_Rational:
            if value != 0:
                return False
        # else a denominator vanished at this point: try another
    numerator, _ = expr.as_numer_denom()
    return symengine.expand(numerator) == 0


class SymEngineExpression(CommutativeRingElement):
    """
    A SymEngine expression as a Sage ring element. Arithmetic is structural (no expansion, no
    cancellation); ``bool``, ``is_zero`` and comparison with ``0`` test for the literal zero.
    ``==`` against anything else is semantic: structural first, then evaluation at random integer
    points (which settles the nonzero case cheaply), and only an expression vanishing at every point is
    expanded.

    EXAMPLES::

        sage: from schubmult.sage.symengine_ring import SymEngineRing
        sage: E = SymEngineRing(ZZ); a = E.variable('a', 1); b = E.variable('b', 1)
        sage: (a + b)^2 - (a + b)^2
        0
        sage: (a + b)^2 - (a^2 + 2*a*b + b^2)
        -(2*a_1*b_1 + a_1**2 + b_1**2) + (a_1 + b_1)**2
        sage: bool(_), _ == 0, _.is_zero(), _.is_identically_zero()
        (True, False, False, True)
        sage: (a + b)^2 == a^2 + 2*a*b + b^2
        True
        sage: hash(a + b) == hash(b + a)
        True
    """

    def __init__(self, parent, expr):
        CommutativeRingElement.__init__(self, parent)
        self._expr = expr

    def expr(self):
        """The underlying SymEngine expression."""
        return self._expr

    def _new(self, expr):
        return type(self)(self.parent(), expr)

    def _repr_(self):
        return str(self._expr)

    def _latex_(self):
        from schubmult.symbolic import latex

        return latex(self._expr._sympy_())

    def __hash__(self):
        return hash(self._expr)

    def __bool__(self):
        return self._expr != 0

    def is_zero(self):
        return not self

    def _richcmp_(self, other, op):
        # against the literal zero (``c != 0`` in printing, ``is_zero``) the test is structural, so that
        # displaying an element never expands anything; otherwise it is semantic
        if op not in (op_EQ, op_NE):
            raise TypeError("SymEngine expressions are not ordered")
        eq = self._expr == other._expr
        if not eq and self._expr != 0 and other._expr != 0:
            eq = is_identically_zero(self._expr - other._expr)
        return eq if op == op_EQ else not eq

    def is_identically_zero(self):
        """Whether the expression is zero as a rational function (``is_zero`` tests for the literal zero)."""
        return is_identically_zero(self._expr)

    def _add_(self, other):
        return self._new(self._expr + other._expr)

    def _sub_(self, other):
        return self._new(self._expr - other._expr)

    def _mul_(self, other):
        return self._new(self._expr * other._expr)

    def _neg_(self):
        return self._new(-self._expr)

    def _div_(self, other):
        # structural quotient (a negative power), as in the double Grothendieck kernels
        if not other:
            raise ZeroDivisionError("division by zero")
        return self._new(self._expr / other._expr)

    def __pow__(self, n, modulus=None):  # noqa: ARG002
        return self._new(self._expr ** int(n))

    def __invert__(self):
        return self.parent().one() / self

    def is_unit(self):
        if self._expr.is_Number:
            return bool(self) and (self.parent().base() is QQ or self._expr in (1, -1))
        return False

    def expand(self):
        """The expanded expression (sum of monomials); exponential in the size of the tree."""
        import symengine

        return self._new(symengine.expand(self._expr))

    def factor(self):
        """The expression in factored form (via SymPy)."""
        import symengine
        import sympy

        return self._new(symengine.sympify(sympy.factor(self._expr._sympy_())))

    def variables(self):
        """The variables occurring in the expression, as ring elements."""
        return tuple(self._new(s) for s in sorted(self._expr.free_symbols, key=str))

    def subs(self, in_dict=None, **kwds):
        """Substitute; keys may be ring elements, SymEngine symbols, or variable names (``y_1``)."""
        import symengine

        d = dict(in_dict or {})
        d.update(kwds)
        P = self.parent()
        return self._new(self._expr.subs({symengine.Symbol(k) if isinstance(k, str) else P(k)._expr: P(v)._expr for k, v in d.items()}))

    def to_polynomial(self, B, named=None):
        """
        The expression as an element of the Sage ring ``B`` (an infinite or finite polynomial ring,
        or the fraction field of one). The 1-based variable ``y_i`` becomes the 0-based ``y_{i-1}``.
        ``named`` maps unindexed symbol names to elements of ``B``.
        """
        return symengine_to_base_ring([self._expr], B, named)[0]

    def _sympy_(self):
        return self._expr._sympy_()


class SymEngineRing(UniqueRepresentation, Parent):
    r"""
    The ring of SymEngine expressions with rational (``scalars=QQ``) or integer (``ZZ``) scalars.

    This is a commutative ring in Sage's sense -- it is the polynomial ring in all the indexed
    variables, presented without a normal form -- with the arithmetic cost model of the kernels:
    multiplication of coefficients is concatenation of factor lists.

    EXAMPLES::

        sage: from schubmult.sage.symengine_ring import SymEngineRing
        sage: E = SymEngineRing(ZZ); E
        Ring of SymEngine expressions over Integer Ring
        sage: E is SymEngineRing(ZZ), E is loads(dumps(E))
        (True, True)
        sage: TestSuite(E).run()
        sage: E.variable('q', 2) * 3 + 1
        1 + 3*q_2
        sage: E(ZZ(1)/2)
        Traceback (most recent call last):
        ...
        TypeError: 1/2 is not an integer
    """

    Element = SymEngineExpression

    def __init__(self, scalars=QQ):
        if scalars is not ZZ and scalars is not QQ:
            raise ValueError("the scalars must be ZZ or QQ (SymEngine numbers)")
        Parent.__init__(self, base=scalars, category=CommutativeRings())

    def _repr_(self):
        return f"Ring of SymEngine expressions over {self.base()}"

    def _element_constructor_(self, x):
        import symengine

        if isinstance(x, SymEngineExpression):
            return self.element_class(self, x._expr)
        if isinstance(x, symengine.Basic):
            return self.element_class(self, x)
        if isinstance(x, int | Integer):
            return self.element_class(self, symengine.Integer(int(x)))
        if isinstance(x, Rational | Fraction):
            num, den = int(x.numerator()), int(x.denominator())
            if den != 1 and self.base() is ZZ:
                raise TypeError(f"{x} is not an integer")
            return self.element_class(self, symengine.Integer(num) if den == 1 else symengine.Rational(num, den))
        if isinstance(x, FractionFieldElement):
            return self(x.numerator()) / self(x.denominator())
        if isinstance(x, Polynomial):
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

            x = PolynomialRing(x.base_ring(), 1, x.parent().variable_name())(x)
        if isinstance(x, MPolynomial | InfinitePolynomial):
            names = (x.polynomial() if isinstance(x, InfinitePolynomial) else x).parent().variable_names()
            named = {n: symengine.Symbol(n) for n in names if parse_sage_name(n) is None}
            return self.element_class(self, sage_polynomial_to_symengine(x, {}, named))
        if isinstance(x, str) or hasattr(x, "_sympy_"):
            return self.element_class(self, symengine.sympify(x))
        raise TypeError(f"cannot make a SymEngine expression from {x!r}")

    def _coerce_map_from_(self, S):
        if S is int or S is ZZ or (S is QQ and self.base() is QQ):
            return True
        if isinstance(S, SymEngineRing):
            return self.base().has_coerce_map_from(S.base())
        if isinstance(S, MPolynomialRing_base | PolynomialRing_generic | InfinitePolynomialRing_sparse):
            return self.base().has_coerce_map_from(S.base_ring())
        if isinstance(S, FractionField_generic):
            return self.has_coerce_map_from(S.ring())
        return None

    def variable(self, letter, i):
        """The schubmult variable ``<letter>_<i>`` (1-based ``i``) as a ring element."""
        from schubmult.symbolic.poly.variables import GeneratingSet

        return self.element_class(self, GeneratingSet(letter)[i])

    def symbol(self, name):
        """An unindexed symbol (e.g. ``'beta'``) as a ring element."""
        import symengine

        return self.element_class(self, symengine.Symbol(name))

    @cached_method
    def zero(self):
        return self(0)

    @cached_method
    def one(self):
        return self(1)

    def _an_element_(self):
        return self.variable("y", 1) + 1

    def some_elements(self):
        y1, y2, y3 = (self.variable("y", i) for i in (1, 2, 3))
        return [self.zero(), self.one(), y1, y1 - y2, (y1 - y3) * (y2 - y3), (y1 + 1) ** 2]

    def is_exact(self):
        return True

    def is_field(self, proof=True):  # noqa: ARG002
        return False

    def is_integral_domain(self, proof=True):  # noqa: ARG002
        return True

    def characteristic(self):
        return ZZ.zero()
