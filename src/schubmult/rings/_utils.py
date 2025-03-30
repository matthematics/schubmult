from functools import cache

from sympy import Basic, Indexed, IndexedBase, sympify


class Zero(Basic):
    def __getitem__(self,i):
        return sympify(0)
    
zero = sympify(0)
zeroindexer = Zero()

NoneVar = 1e10
ZeroVar = 0

class ZeroVarClass(IndexedBase):
    base = ZeroVar
    args = (0,0)

    def __new__(cls):
        obj = IndexedBase.__new__(cls,zeroindexer)
        return obj

    def __getitem__(self,i):
        return zero


class NoneVarClass(IndexedBase):
    base = NoneVar
    args = (0,0)

    def __new__(cls):
        obj = IndexedBase.__new__(cls,zeroindexer)
        return zero

    def __getitem__(self,i):
        return 0



@cache
def poly_ring(v: str):
    if v == ZeroVar:
        return tuple([0 for i in range(100)])
    if v == NoneVar:
        return tuple([0 for i in range(100)])
    return IndexedBase(v)
