from functools import cache
from symengine import symarray

@cache
def poly_ring(v: str):
    return tuple(symarray(v, 100).tolist())