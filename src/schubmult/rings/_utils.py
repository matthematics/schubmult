from functools import cache

from symengine import symarray
from sympy import Dummy

NoneVar = 0 #Dummy("None")
ZeroVar = 0

@cache
def poly_ring(v: str):
    if v == ZeroVar or v == NoneVar:
        return tuple([0 for i in range(100)])    
    return tuple(symarray(v, 100).tolist())
