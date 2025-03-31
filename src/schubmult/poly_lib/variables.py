# class generators with base
# symbols cls argument!

import sympy
from symengine import symbols
from symengine.lib.symengine_wrapper import Symbol


class ISymbol(Symbol):
    def __new__(cls,name):
        obj = Symbol.__new__(cls,name)
        obj.base =None
        obj.index = 0
        return obj
    
    def __hash__(self):
        return hash(self.name)
    
    def __eq__(self, other):
        if not isinstance(other, Symbol):
            return False
        return other.name == self.name

class GeneratingSet:
    def __init__(self, name):
        self._symbols_arr = symbols(f"{name}(0:100)",cls=ISymbol)
        self.label = name
        for i in range(len(self._symbols_arr)):
            self._symbols_arr[i].base = self
            self._symbols_arr[i].index = i

    def __getitem__(self, i):
        return self._symbols_arr[i]

    def __hash__(self):
        return hash(self._symbols_arr)

    def __eq__(self,other):
        return isinstance(other,GeneratingSet) and self.label==other.label

def is_indexed(x):
    return hasattr(x, "index")

