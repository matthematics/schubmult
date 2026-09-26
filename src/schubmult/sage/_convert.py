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
