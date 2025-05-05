import symengine
import sympy

from .symmetric_polynomials.complete_sym import CompleteSym
from .symmetric_polynomials.elem_sym import ElemSym


def sympify(val):
    try:
        return symengine.sympify(val)
    except symengine.SympifyError:
        # print(f"Bagels {val=}")
        return sympy.sympify(val)


def expand(val, **kwargs):
    try:
        return symengine.expand(val, **kwargs)
    except Exception:
        # print(f"Bogels {val=}")
        return sympy.expand(val, **kwargs)
