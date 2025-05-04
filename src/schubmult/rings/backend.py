import symengine
import sympy

from .symmetric_polynomials.symengine.elem_sym import ElemSym
from .symmetric_polynomials.sympy.complete_sym import CompleteSym


def sympify(val):
    try:
        return symengine.sympify(val)
    except symengine.SympifyError:
        print(f"Bagels {val=}")
        return sympy.sympify(val)


def expand(val, **kwargs):
    try:
        return symengine.expand(val, **kwargs)
    except Exception:
        print(f"Bogels {val=}")
        return sympy.expand(val, **kwargs)
