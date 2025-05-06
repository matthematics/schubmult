import symengine.lib.symengine_wrapper as sw
from sympy.core.expr import Expr

class Test(sw.Symbol, Printable):

    def __init__(self, expr, *args):
        self._base_args = args
        super().__init__(self, store_pickle=True)

    def encode(self, *args):
        from sympy.printing.str import sstr
        return sstr(self).encode(*args)

    @property
    def args(self):
        return self._base_args
    
    def __reduce__(self):
        return (self.__class__, self.args)

    def _sympystr(self, printer):
        return printer._print(self.args)
    
    def _sympyrepr(self, printer):
        return printer._print("Felshnoop")
    
from schubmult import *

a = Test("a",12,43)
b = Test(13,1)
print(a*b)
print(type((a*b).args[1]))