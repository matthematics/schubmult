def sstr(expr, **settings):
    p = SchubStrPrinter(settings)
    s = p.doprint(expr)
    return s

import sympy
from sympy.printing.str import StrPrinter


class SchubStrPrinter(StrPrinter):

    def _print_Mul(self, expr):
        print("bagels")
        expr = sympy.Mul(*expr.args)
        return super()._print_Mul(expr)

