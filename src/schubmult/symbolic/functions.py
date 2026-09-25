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
