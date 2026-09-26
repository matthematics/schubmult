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


def sage_polynomial_to_symengine(p, gensets, named=None):
    """Sage polynomial (finite or infinite polynomial ring) -> SymEngine expression.

    ``gensets`` maps a letter to a schubmult ``GeneratingSet``; the Sage variable ``a<i>``/``a_<i>``
    becomes ``gensets[a][i + 1]``. Letters not in ``gensets`` are created on the fly. ``named`` maps
    unindexed Sage variable names (``'beta'``) to SymEngine symbols.
    """
    import symengine
    from sage.rings.polynomial.infinite_polynomial_element import InfinitePolynomial

    from schubmult.symbolic.poly.variables import GeneratingSet

    named = named or {}
    if isinstance(p, InfinitePolynomial):
        p = p.polynomial()
    names = p.parent().variable_names()
    syms = []
    for name in names:
        if name in named:
            syms.append(named[name])
            continue
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


def symengine_to_base_ring(exprs, B, named=None):
    """Batch-convert SymEngine expressions to elements of the coefficient ring ``B``.

    ``B`` may be an ``InfinitePolynomialRing``, the fraction field of one (double Grothendieck
    coefficients are rational in ``beta`` and ``y``), or a finite ``PolynomialRing`` (e.g. ``R[beta]``).
    ``named`` maps unindexed SymEngine symbol names to elements of the finite ring underneath ``B``
    (``{'\u03b2': beta}``).

    A generic tree walk doing every ``+``/``*`` through ``InfinitePolynomial`` arithmetic is Python-level
    and ~10x slower than the libsingular ring underneath; here the trees are evaluated directly in the
    underlying finite ring (grown first to cover every index that occurs) and wrapped without
    conversion. Coefficients are never symbolically expanded on the schubmult side; the normal form is
    computed by libsingular. Coefficients of one product share most of their subtrees (the same
    ``(y_i - z_j)`` factors and partial products recur), so values are memoized on the SymEngine node
    across the whole batch.
    """
    import symengine
    from sage.rings.fraction_field import FractionField_generic
    from sage.rings.polynomial.infinite_polynomial_ring import InfinitePolynomialRing_sparse
    from sage.rings.rational_field import QQ

    named = named or {}
    exprs = [symengine.sympify(e) for e in exprs]
    fraction_field = B if isinstance(B, FractionField_generic) else None
    D = fraction_field.ring() if fraction_field is not None else B

    if isinstance(D, InfinitePolynomialRing_sparse):  # dense is a subclass
        letters = {name: i for i, name in enumerate(D.variable_names())}
        max_index = dict.fromkeys(letters, 0)
        for e in exprs:
            for s in e.free_symbols:
                m = _SCHUB_NAME.match(str(s))
                if m is not None and m.group(1) in letters:
                    max_index[m.group(1)] = max(max_index[m.group(1)], int(m.group(2)) - 1)
                elif str(s) not in named:
                    raise ValueError(f"cannot convert symbol {s} to an element of {B}")
        for letter, idx in max_index.items():
            D.gen(letters[letter])[idx]
        P = D.gen(0)[0].polynomial().parent()  # the (now large enough) libsingular ring
        poly_class = type(D.gen(0)[0])

        def wrap(p):
            return poly_class(D, p)
    else:
        P = D

        def wrap(p):
            return p

    gens = dict(zip(P.variable_names(), P.gens()))
    named_in_P = {name: P(v) for name, v in named.items()}
    memo = {}
    zero, one = P.zero(), P.one()

    # Values are polynomials in P, or flat fractions ``(numerator, {atom: exponent})`` once a negative
    # power has been met, the denominator being ``prod atom**exponent``.  This mirrors the kernels
    # (``grothmult_double`` keeps its coefficients over ``prod (1 + beta*y_i)**e``): sums go over the
    # lcm of the atom dicts, so intermediate results never need a gcd -- the finite fraction field
    # would reduce at every ``+``/``*``, and those gcds of the big unexpanded numerators dominate
    # everything else.  One gcd per coefficient at the end reduces the result.
    def mul(a, b):
        if type(a) is tuple or type(b) is tuple:
            an, ad = a if type(a) is tuple else (a, {})
            bn, bd = b if type(b) is tuple else (b, {})
            if not bd:
                return (an * bn, ad)
            if not ad:
                return (an * bn, bd)
            d = dict(ad)
            for atom, e in bd.items():
                d[atom] = d.get(atom, 0) + e
            return (an * bn, d)
        return a * b

    def add(a, b):
        if type(a) is tuple or type(b) is tuple:
            an, ad = a if type(a) is tuple else (a, {})
            bn, bd = b if type(b) is tuple else (b, {})
            if ad == bd:
                return (an + bn, ad)
            d = dict(ad)
            for atom, e in bd.items():
                if e > d.get(atom, 0):
                    d[atom] = e
            for atom, e in d.items():
                ea, eb = ad.get(atom, 0), bd.get(atom, 0)
                if ea < e:
                    an = an * atom ** (e - ea)
                if eb < e:
                    bn = bn * atom ** (e - eb)
            return (an + bn, d)
        return a + b

    def power(a, k):
        if k < 0:
            an, ad = a if type(a) is tuple else (a, {})
            n = one
            for atom, e in ad.items():
                n = n * atom ** (e * -k)
            return (n, {an: -k})
        if type(a) is tuple:
            return (a[0] ** k, {atom: e * k for atom, e in a[1].items()})
        return a**k

    def reduce_fraction(v):
        """``(numerator, denominator)`` in P, reduced (one gcd, over the scalars: cheap)."""
        if type(v) is not tuple:
            return v, one
        n, d = v
        den = one
        for atom, e in d.items():
            den = den * atom**e
        g = n.gcd(den)
        if not g.is_one():
            n, den = n // g, den // g
        return n, den

    def go(e):
        v = memo.get(e)
        if v is not None:
            return v
        if e.is_Integer:
            v = P(int(e))
        elif e.is_Symbol:
            name = str(e)
            if name in named_in_P:
                v = named_in_P[name]
            else:
                m = _SCHUB_NAME.match(name)
                v = gens[f"{m.group(1)}_{int(m.group(2)) - 1}"]
        elif e.is_Add:
            v = zero
            for a in e.args:
                v = add(v, go(a))
        elif e.is_Mul:
            v = one
            for a in e.args:
                v = mul(v, go(a))
        elif e.is_Pow:
            base, exp = e.args
            if not exp.is_Integer or (int(exp) < 0 and fraction_field is None):
                raise ValueError(f"cannot convert {e} to an element of {B}")
            v = power(go(base), int(exp))
        elif e.is_Rational:
            v = P(QQ((int(e.p), int(e.q))))
        else:
            raise ValueError(f"cannot convert {e} ({type(e).__name__}) to Sage")
        memo[e] = v
        return v

    if fraction_field is None:
        return [wrap(go(e)) for e in exprs]
    frac_class = fraction_field._element_class
    out = []
    for e in exprs:
        n, d = reduce_fraction(go(e))
        out.append(frac_class(fraction_field, wrap(n), wrap(d), coerce=False, reduce=False))
    return out


def symengine_to_sage(expr, variable, scalar, named=None):
    """SymEngine expression -> Sage element.

    ``variable(letter, i)`` returns the Sage element for the schubmult symbol ``letter_i`` (1-based ``i``);
    ``scalar(n)`` converts a Python ``int``/``Fraction``-like rational to the target ring; ``named`` maps
    unindexed symbol names (``'\u03b2'``) to target elements.
    """
    from fractions import Fraction

    named = named or {}

    def go(e):
        if e.is_Integer:
            return scalar(int(e))
        if e.is_Rational:
            return scalar(Fraction(int(e.p), int(e.q)))
        if e.is_Symbol:
            if str(e) in named:
                return named[str(e)]
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
