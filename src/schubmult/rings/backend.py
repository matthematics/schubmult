import symengine
import sympy


def sympify(val):
    try:
        return symengine.sympify(val)
    except symengine.SympifyError:
        return sympy.sympify(val)


def expand(val, **kwargs):
    try:
        return symengine.expand(val, **kwargs)
    except Exception:
        return sympy.expand(val, **kwargs)
