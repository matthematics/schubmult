import symengine
import sympy

from .symbol import SymengineSymbol


def expand(obj, **kwargs):
    if len(kwargs.keys()):
        return sympy.expand(obj, **kwargs)
    return symengine.expand(obj)

S = symengine.S

def symbols(*args, **kwargs):
    return symengine.symbols(*args, **kwargs)