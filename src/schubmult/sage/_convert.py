"""Conversions between schubmult's SymEngine expressions and Sage ring elements.

Index conventions differ by one: schubmult's generating sets are 1-indexed (``x_1`` is the first
variable, ``x_0`` unused) while Sage's Schubert polynomials expand into ``x0, x1, ...``. Every
conversion here shifts accordingly, so ``y_3`` on the schubmult side is ``y_2`` / ``y2`` on the Sage side.
"""

import re

_SAGE_NAME = re.compile(r"^([A-Za-z]+)_?([0-9]+)$")
_SCHUB_NAME = re.compile(r"^([^_]+)_([0-9]+)$")


def parse_sage_name(name):
    """``'x3'`` or ``'y_3'`` -> ``('x', 3)`` / ``('y', 3)`` (0-based Sage index); ``None`` if not indexed."""
    m = _SAGE_NAME.match(str(name))
    return (m.group(1), int(m.group(2))) if m else None


def sage_coefficient_to_symengine(c):
    """Base-ring scalar (integer or rational) to a SymEngine number."""
    import symengine

    if hasattr(c, "numerator") and hasattr(c, "denominator"):
        num, den = int(c.numerator()), int(c.denominator())
        return symengine.Integer(num) if den == 1 else symengine.Rational(num, den)
    return symengine.sympify(str(c))


def sage_polynomial_to_symengine(p, gensets):
    """Sage polynomial (finite or infinite polynomial ring) -> SymEngine expression.

    ``gensets`` maps a letter to a schubmult ``GeneratingSet``; the Sage variable ``a<i>``/``a_<i>``
    becomes ``gensets[a][i + 1]``. Letters not in ``gensets`` are created on the fly.
    """
    import symengine
    from sage.rings.polynomial.infinite_polynomial_element import InfinitePolynomial

    from schubmult.symbolic.poly.variables import GeneratingSet

    if isinstance(p, InfinitePolynomial):
        p = p.polynomial()
    names = p.parent().variable_names()
    syms = []
    for name in names:
        parsed = parse_sage_name(name)
        if parsed is None:
            raise ValueError(f"variable {name!r} is not of the form <letter><index>")
        letter, idx = parsed
        gs = gensets.get(letter)
        if gs is None:
            gs = gensets[letter] = GeneratingSet(letter)
        syms.append(gs[idx + 1])
    result = symengine.Integer(0)
    for exps, c in p.dict().items():
        term = sage_coefficient_to_symengine(c)
        for s, e in zip(syms, exps):
            if e:
                term *= s**e
        result += term
    return result


def symengine_to_infinite_polynomials(exprs, B):
    """Batch-convert SymEngine expressions to elements of the ``InfinitePolynomialRing`` ``B``.

    The generic tree walk does every ``+``/``*`` through ``InfinitePolynomial`` arithmetic, which is
    Python-level and ~10x slower than the libsingular ring underneath. Here the trees are evaluated
    directly in ``B``'s underlying finite ring (grown first to cover every index that occurs) and the
    results wrapped without conversion. Coefficients are never symbolically expanded on the schubmult
    side; the polynomial normal form is computed by libsingular.
    """
    import symengine
    from sage.rings.rational_field import QQ

    exprs = [symengine.sympify(e) for e in exprs]
    letters = {name: i for i, name in enumerate(B.variable_names())}
    # grow the underlying ring so it has every variable we will need (Sage index = schubmult index - 1)
    max_index = dict.fromkeys(letters, 0)
    for e in exprs:
        for s in e.free_symbols:
            m = _SCHUB_NAME.match(str(s))
            if m is None or m.group(1) not in letters:
                raise ValueError(f"cannot convert symbol {s} to an element of {B}")
            max_index[m.group(1)] = max(max_index[m.group(1)], int(m.group(2)) - 1)
    for letter, idx in max_index.items():
        B.gen(letters[letter])[idx]
    P = B.gen(0)[0].polynomial().parent()  # the (now large enough) libsingular ring
    gens = dict(zip(P.variable_names(), P.gens()))
    element_class = type(B.gen(0)[0])
    # Coefficients of one product share most of their subtrees (the same (y_i - z_j) factors and
    # partial products recur across terms), so memoizing on the SymEngine node cuts the walk ~10x.
    memo = {}
    zero, one = P.zero(), P.one()

    def go(e):
        v = memo.get(e)
        if v is not None:
            return v
        if e.is_Integer:
            v = P(int(e))
        elif e.is_Symbol:
            m = _SCHUB_NAME.match(str(e))
            v = gens[f"{m.group(1)}_{int(m.group(2)) - 1}"]
        elif e.is_Add:
            v = zero
            for a in e.args:
                v += go(a)
        elif e.is_Mul:
            v = one
            for a in e.args:
                v *= go(a)
        elif e.is_Pow:
            base, exp = e.args
            if not exp.is_Integer or int(exp) < 0:
                raise ValueError(f"cannot convert {e} to a polynomial")
            v = go(base) ** int(exp)
        elif e.is_Rational:
            v = P(QQ((int(e.p), int(e.q))))
        else:
            raise ValueError(f"cannot convert {e} ({type(e).__name__}) to Sage")
        memo[e] = v
        return v

    return [element_class(B, go(e)) for e in exprs]


def symengine_to_sage(expr, variable, scalar):
    """SymEngine expression -> Sage element.

    ``variable(letter, i)`` returns the Sage element for the schubmult symbol ``letter_i`` (1-based ``i``);
    ``scalar(n)`` converts a Python ``int``/``Fraction``-like rational to the target ring.
    """
    from fractions import Fraction

    def go(e):
        if e.is_Integer:
            return scalar(int(e))
        if e.is_Rational:
            return scalar(Fraction(int(e.p), int(e.q)))
        if e.is_Symbol:
            m = _SCHUB_NAME.match(str(e))
            if m is None:
                raise ValueError(f"cannot convert symbol {e} to Sage")
            return variable(m.group(1), int(m.group(2)))
        if e.is_Add:
            out = scalar(0)
            for a in e.args:
                out += go(a)
            return out
        if e.is_Mul:
            out = scalar(1)
            for a in e.args:
                out *= go(a)
            return out
        if e.is_Pow:
            base, exp = e.args
            if not exp.is_Integer or int(exp) < 0:
                raise ValueError(f"cannot convert {e} to a polynomial")
            return go(base) ** int(exp)
        raise ValueError(f"cannot convert {e} ({type(e).__name__}) to Sage")

    import symengine

    return go(symengine.sympify(expr))
