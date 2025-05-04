from sympy import Wild, sympify

from schubmult.utils.ring_utils import NotEnoughGeneratorsError

from .elem_sym import ElemSym


def split_out_vars(expr, vars1, vars2):
    expr = sympify(expr)
    if not expr.args:
        return expr
    if hasattr(expr, "split_out_vars"):
        try:
            return expr.split_out_vars(vars1, vars2)
        except NotEnoughGeneratorsError:
            return expr
    return expr.func(*[split_out_vars(arg, vars1, vars2) for arg in expr.args])


def pull_out_vars(expr, var1, var2, min_degree=1):
    expr = sympify(expr)
    if hasattr(expr, "pull_out_vars"):
        return expr.pull_out_vars(var1, var2, min_degree)
    if not expr.args:
        return expr
    return expr.func(*[pull_out_vars(arg, var1, var2, min_degree) for arg in expr.args])

def elem_sym_unify(expr, arg=None):
    expr = sympify(expr)

    if not expr.args:
        return expr
    if isinstance(expr, ElemSym):
        return expr

    if not arg:
        for arg in expr.args:
            expr = elem_sym_unify(expr, arg)
        return expr.func(*[elem_sym_unify(arg) for arg in expr.args])
    expr2 = expr
    if isinstance(arg, ElemSym):
        v1 = Wild("v_1")
        v2 = Wild("v_2")
        rep_pattern = ElemSym(arg._p, arg._k + 1, [*arg.genvars, v1], [*arg.coeffvars, v2])
        pattern = arg + ElemSym(1, 1, [v1], [v2]) * ElemSym(arg._p - 1, arg._k, arg.genvars, [*arg.coeffvars, v2])
        expr2 = expr.replace(pattern, rep_pattern)
    if arg.args:
        for arg2 in arg.args:
            expr2 = elem_sym_unify(expr2, arg2)
    if expr != expr2:
        return elem_sym_unify(expr2)
    return expr
