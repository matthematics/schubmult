from functools import cache

import symengine
from sympy import Basic, Expr, Integer, Mul, sympify

import schubmult.poly_lib.variables as vsyme


# Have to create add and mul handlers?
class ISymbol(Expr):
    
    _op_priority = 124135135235
    def __new__(cls, base, index):
        return ISymbol.__xnew_cached__(cls, base, index)

    is_commutative = True

    @staticmethod
    @cache
    def __xnew_cached__(_class, base, index):
        return ISymbol.__xnew__(_class, base, index)

    @staticmethod
    def __xnew__(_class, base, index):
        return Expr.__new__(_class, base, Integer(index))

    def _sympystr(self, printer):
        return printer._print_Symbol(self)

    def _latex(self, printer):
        return printer._print_Symbol(self)

    # def __mul__(self, other):
    #     return SMul(self, other)
    
    # def __rmul__(self, other):
    #     return SMul(other, self)

    # @property
    # def _mul_handler(self):
    #     return SMul

    @property
    def base(self):
        return self.args[0]

    @property
    def name(self):
        return f"{self.args[0].label}_{self.args[1]}"

    @property
    def index(self):
        return self.args[1]

    @property
    def indices(self):
        return (self.args[1],)

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        try:
            return self.name == sympify(other).name
        except Exception:
            return NotImplemented
 
    def _symengine_(self):
            return vsyme.ISymbol(self.name,base=self.args[0],index=self.args[1])
    
    
class SMul(Mul,ISymbol):
    
    def __new__(cls, *args):
        obj = Mul.__new__(cls,*args)
        return obj
    
    def _symengine_(self):
        return vsyme.SMul(*[symengine.sympify(arg) for arg in self.args])
    
Basic._constructor_postprocessor_mapping[ISymbol] = {Mul: lambda expr: SMul(*expr.args)}