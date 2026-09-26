"""SymEngine-first wrappers (`expand`, `symbols`, `sympify`) that fall back to SymPy, plus small helpers."""

import operator
from functools import reduce

import symengine
import symengine.lib.symengine_wrapper as sw


def expand(obj, **kwargs):
    """Expand with SymEngine; use SymPy if keyword options are given or SymEngine fails."""
    if len(kwargs.keys()):
        import sympy

        return symengine.sympify(sympy.expand(obj, **kwargs))
    try:
        return symengine.expand(obj)
    except Exception:
        import sympy

        return sympy.expand(obj)


def symbols(*args, **kwargs):
    """SymEngine ``symbols``."""
    return symengine.symbols(*args, **kwargs)


def sympify(val):
    """SymEngine ``sympify``, falling back to SymPy for objects SymEngine cannot convert."""
    try:
        return symengine.sympify(val)
    except symengine.SympifyError:
        import sympy

        return sympy.sympify(val)


def is_of_func_type(elem, typ):
    """``isinstance`` that also sees through SymEngine ``PyFunction`` wrappers around SymPy functions."""
    return isinstance(elem, typ) or (isinstance(elem, sw.PyFunction) and isinstance(elem.pyobject(), typ))


def expand_seq(seq, genset):
    """The monomial ``genset[1]**seq[0] * genset[2]**seq[1] * ...`` (1-indexed generators)."""
    return prod([genset[i + 1] ** seq[i] for i in range(len(seq))])


def prod(a, start=1):
    """Product of the elements of ``a`` times ``start`` (same as ``sympy.prod``)."""
    return reduce(operator.mul, a, start)


def efficient_subs(expr, subs_dict):
    """``expr.subs`` restricted to the entries of ``subs_dict`` that actually occur in ``expr``."""
    expr = sympify(expr)
    subs_dict_new = {s: subs_dict[s] for s in expr.free_symbols if s in subs_dict}
    return expr.subs(subs_dict_new)


def vanish_at_random_points(exprs, trials=2, seed=1, bound=10**6):
    """For each expression, whether it evaluates to zero at ``trials`` random integer points (exact arithmetic).

    A polynomial that is identically zero always does; a nonzero one of modest degree vanishes at a
    random point of ``[-bound, bound]^n`` with negligible probability, so this is a fast surrogate for
    ``expand(e) == 0`` that never touches the (possibly enormous) expanded form. Values of shared
    subtrees are memoized across the whole batch, which is what makes it cheap: the coefficients of one
    Schubert product reuse the same ``(y_i - z_j)`` factors and partial products over and over.

    Expressions containing nodes other than numbers, symbols, ``Add``, ``Mul`` and integer ``Pow`` fall
    back to ``expand(e) == 0``.
    """
    import random
    from fractions import Fraction

    rng = random.Random(seed)
    values = {}
    memo = {}

    def go(e):
        r = memo.get(e)
        if r is not None:
            return r
        if e.is_Integer:
            r = (int(e),) * trials
        elif e.is_Rational:
            r = (Fraction(int(e.p), int(e.q)),) * trials
        elif e.is_Symbol:
            r = values.get(e)
            if r is None:
                r = values[e] = tuple(rng.randint(-bound, bound) for _ in range(trials))
        elif e.is_Add:
            acc = [0] * trials
            for a in e.args:
                for i, x in enumerate(go(a)):
                    acc[i] += x
            r = tuple(acc)
        elif e.is_Mul:
            acc = [1] * trials
            for a in e.args:
                for i, x in enumerate(go(a)):
                    acc[i] *= x
            r = tuple(acc)
        elif e.is_Pow and e.args[1].is_Integer:
            k = int(e.args[1])
            r = tuple(Fraction(x) ** k if k < 0 else x**k for x in go(e.args[0]))
        else:
            raise TypeError(f"unsupported node {type(e).__name__}")
        memo[e] = r
        return r

    out = []
    for e in exprs:
        try:
            out.append(not any(go(symengine.sympify(e))))
        except (TypeError, symengine.SympifyError):
            out.append(expand(e) == 0)
    return out
