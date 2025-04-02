from functools import cache

from symengine import sympify

from schubmult.poly_lib import GeneratingSet

NoneVar = 1e10
ZeroVar = 0





@cache
def poly_ring(v: str):
    if v == ZeroVar:
        return tuple([sympify(0) for i in range(100)])
    if v == NoneVar:
        return tuple([sympify(0) for i in range(100)])
    return GeneratingSet(str(v))
