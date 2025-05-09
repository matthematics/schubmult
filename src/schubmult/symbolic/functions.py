import symengine
import sympy

from .symbol import SymengineSymbol


def expand(obj, **kwargs):
    if len(kwargs.keys()):
        return sympy.expand(obj, **kwargs)
    return symengine.expand(obj)

S = symengine.S

def symbols(*args, **kwargs):
    if "cls" not in kwargs:
        return sympy.symbols(*args, cls=SymengineSymbol, **kwargs)
    return sympy.symbols(*args, **kwargs)