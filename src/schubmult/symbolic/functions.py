import symengine
import sympy

from .base_printing import sstr

# Add = symengine.Add
# Mul = symengine.Mul
# Pow = symengine.Pow
# SympifyError = symengine.SympifyError
# Symbol = symengine.Symbol
# S = symengine.S


def expand(obj, **kwargs):
    if len(kwargs.keys()):
        return symengine.sympify(sympy.expand(obj, **kwargs))
    return symengine.expand(obj)
    # except Exception:
    #     return symengine.sympify(sympy.expand(obj))


def symbols(*args, **kwargs):
    return symengine.symbols(*args, **kwargs)


def sympify(val):
    # try:
    return symengine.sympify(val)
    # except symengine.SympifyError:
    #     # print(f"Bagels {val=}")
    #     return sympy.sympify(val)
