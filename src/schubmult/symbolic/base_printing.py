import sympy
from sympy.printing.str import StrPrinter


class SchubStrPrinter(StrPrinter):
    def _print_Add(self, expr):
        try:
            expr = sympy.sympify(sympy.Add(*[sympy.sympify(arg) for arg in expr.args]))
        except Exception:
            pass
        return super()._print_Add(sympy.sympify(expr))

    def _print_Mul(self, expr):
        expr = sympy.Mul(*[sympy.sympify(arg) for arg in expr.args])
        return super()._print_Mul(expr)

    # def _print(self, expr):
    #     return super()._print(sympy.sympify(expr))


def sstr(expr, **settings):
    p = SchubStrPrinter(settings)
    s = p.doprint(expr)
    return s
