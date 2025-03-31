from functools import cache

from sympy import IndexedBase

NoneVar = 1e10
ZeroVar = 0





@cache
def poly_ring(v: str):
    if v == ZeroVar:
        return tuple([0 for i in range(100)])
    if v == NoneVar:
        return tuple([0 for i in range(100)])
    return IndexedBase(v)
