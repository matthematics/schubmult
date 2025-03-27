from functools import cache

from symengine import symarray


@cache
def poly_ring(v: str):
    if v is None:
        return tuple([0 for i in range(100)])
    return tuple(symarray(v, 100).tolist())
