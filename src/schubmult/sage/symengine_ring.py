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

from collections import defaultdict
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

_PRIME = (1 << 61) - 1
_TRIAL_POINTS = 3
REPR_TREE_LIMIT = 400  # larger coefficient trees print in normal form
_random = Random(0x5EED)


def _evaluator_mod_p():
    """A memoized evaluator of SymEngine expressions at a random point modulo ``_PRIME``.

    Returns ``None`` for an expression whose denominator vanishes at the point. The memo is keyed by
    node, so the shared subtrees of a batch of coefficients are evaluated once.
    """
    point = defaultdict(lambda: _random.randrange(1, _PRIME))
    memo = {}
    sentinel = object()

    def ev(e):
        v = memo.get(e, sentinel)
        if v is not sentinel:
            return v
        if e.is_Integer:
            v = int(e) % _PRIME
        elif e.is_Symbol:
            v = point[e]
        elif e.is_Add:
            v = 0
            for a in e.args:
                x = ev(a)
                if x is None:
                    v = None
                    break
                v += x
            else:
                v %= _PRIME
        elif e.is_Mul:
            v = 1
            for a in e.args:
                x = ev(a)
                if x is None:
                    v = None
                    break
                v = v * x % _PRIME
        elif e.is_Pow:
            base, exp = e.args
            b, k = ev(base), int(exp)
            v = None if b is None or (k < 0 and b == 0) else pow(b, k, _PRIME)
        elif e.is_Rational:
            v = int(e.p) * pow(int(e.q), -1, _PRIME) % _PRIME
        else:
            raise TypeError(f"cannot evaluate {e} ({type(e).__name__})")
        memo[e] = v
        return v

    return ev


def nonzero_mask(exprs, points=2):
    """For each SymEngine expression, whether it is nonzero at one of ``points`` random points modulo
    the prime `2^{61} - 1` (so certainly nonzero as a rational function).

    One shared evaluation per point over all the expressions: the coefficients of one product share
    most of their subtrees. An expression vanishing at every point is reported zero -- for a nonzero
    rational function of degree `d` that happens with probability at most `(d / 2^{61})^{\text{points}}`;
    one undefined at every point (a denominator vanishing there) is reported nonzero.
    """
    exprs = list(exprs)
    mask = [False] * len(exprs)
    defined = [False] * len(exprs)
    undecided = range(len(exprs))
    for _ in range(points):
        ev = _evaluator_mod_p()
        for i in undecided:
            v = ev(exprs[i])
            if v is not None:
                defined[i] = True
                if v != 0:
                    mask[i] = True
        undecided = [i for i in undecided if not mask[i]]
        if not undecided:
            break
    return [m or not d for m, d in zip(mask, defined)]


def is_identically_zero(expr):
    """Whether the SymEngine expression ``expr`` is zero as a rational function.

    Evaluation at a few random points first (one pass over the shared subtrees): a nonzero value
    proves nonvanishing, which is the common case. Only an expression vanishing at every point is
    expanded (over a common denominator) for a proof.
    """
    import symengine

    if expr == 0:
        return True
    if nonzero_mask([expr], _TRIAL_POINTS)[0]:
        return False
    numerator, _ = expr.as_numer_denom()
    return symengine.expand(numerator) == 0


def symengine_to_sage_scalar(e, R):
    """SymEngine integer or rational -> element of ``R``."""
    return R(int(e)) if e.is_Integer else R(QQ((int(e.p), int(e.q))))


def tree_size(expr, limit):
    """Number of nodes of the expression *tree* (shared subtrees counted each time), stopping past ``limit``."""
    memo = {}

    def go(e):
        v = memo.get(e)
        if v is None:
            v = 1
            for a in e.args:
                v += go(a)
                if v > limit:
                    break
            memo[e] = v
        return v

    return go(expr)


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
        # the kernels' coefficients are DAGs with heavy sharing whose trees run to megabytes of text
        if tree_size(self._expr, REPR_TREE_LIMIT) <= REPR_TREE_LIMIT:
            return str(self._expr)
        return str(self.normal_form())

    def _latex_(self):
        from schubmult.symbolic import latex

        if tree_size(self._expr, REPR_TREE_LIMIT) <= REPR_TREE_LIMIT:
            return latex(self._expr._sympy_())
        return self.normal_form()._latex_()

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

    def __pow__(self, n, modulus=None):
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

    def normal_form(self):
        """
        The expression as a Sage polynomial (or fraction) in a ring on the symbols that occur, named as
        in the expression (``y_3`` stays ``y_3``; ``β`` becomes ``beta``). Computed in libsingular with
        the subtrees shared, which is far cheaper than ``expand``; this is how large coefficients print.

        EXAMPLES::

            sage: from schubmult.sage.symengine_ring import SymEngineRing
            sage: E = SymEngineRing(ZZ); y = [E.variable('y', i) for i in range(4)]
            sage: ((y[1] - y[3]) * (y[2] - y[3])).normal_form()
            y_1*y_2 - y_1*y_3 - y_2*y_3 + y_3^2
            sage: (y[1] / (1 + E.symbol('β') * y[2])).normal_form()
            y_1/(y_2*beta + 1)
        """
        return normal_forms([self._expr], self.parent().base())[0]

    def _sympy_(self):
        return self._expr._sympy_()


def normal_forms(exprs, R):
    """Sage polynomials (or fractions) over ``R`` for a batch of SymEngine expressions, in one ring on
    the union of their symbols and with one memo: the coefficients of a product share their subtrees."""
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

    exprs = list(exprs)
    symbols = sorted({s for e in exprs for s in e.free_symbols}, key=str)
    if not symbols:
        return [symengine_to_sage_scalar(e, R) for e in exprs]
    P = PolynomialRing(R, len(symbols), ["beta" if str(s) == "\u03b2" else str(s) for s in symbols])
    gens = dict(zip(symbols, P.gens()))
    memo = {}
    zero, one = P.zero(), P.one()

    def go(e):
        v = memo.get(e)
        if v is not None:
            return v
        if e.is_Integer or e.is_Rational:
            v = symengine_to_sage_scalar(e, P)
        elif e.is_Symbol:
            v = gens[e]
        elif e.is_Add:
            v = zero
            for a in e.args:
                v = v + go(a)
        elif e.is_Mul:
            v = one
            for a in e.args:
                v = v * go(a)
        elif e.is_Pow:
            base, exp = e.args
            v = go(base) ** int(exp)  # negative: into the fraction field
        else:
            raise TypeError(f"cannot convert {e} ({type(e).__name__}) to Sage")
        memo[e] = v
        return v

    out = []
    for e in exprs:
        v = go(e)
        if isinstance(v, FractionFieldElement) and v.denominator().is_one():
            v = v.numerator()
        out.append(v)
    return out


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
