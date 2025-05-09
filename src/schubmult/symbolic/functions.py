import symengine
import sympy

from .base_printing import sstr
from .symbol import SymengineSymbol

Add = symengine.Add
Mul = symengine.Mul
Pow = symengine.Pow
SympifyError = symengine.SympifyError
Symbol = symengine.Symbol

sympy.init_printing(str_printer = sstr)

def expand(obj, **kwargs):
    if len(kwargs.keys()):
        return sympy.expand(obj, **kwargs)
    #print(f"bagels {type(obj)=} {[(type(arg),type(v)) for arg,v in obj.items()]=}")
    try:
        return symengine.expand(obj)
    except Exception:
        return symengine.sympify(sympy.expand(obj))

S = symengine.S

def symbols(*args, **kwargs):
    return symengine.symbols(*args, **kwargs)


def sympify(val):
    try:
        return symengine.sympify(val)
    except symengine.SympifyError:
        # print(f"Bagels {val=}")
        return sympy.sympify(val)