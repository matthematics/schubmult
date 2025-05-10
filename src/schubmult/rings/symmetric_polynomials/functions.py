from schubmult.rings.variables import NotEnoughGeneratorsError
from schubmult.symbolic import Add, Mul, Pow, expand, is_of_func_type, sympify

from .elem_sym import FactorialElemSym


def genvars(obj):
    try:
        return obj.genvars
    except AttributeError:
        return obj.pyobject().genvars


def coeffvars(obj):
    try:
        return obj.coeffvars
    except AttributeError:
        return obj.pyobject().coeffvars


def degree(obj):
    try:
        return obj._p
    except AttributeError:
        return obj.pyobject()._p


def numvars(obj):
    try:
        return obj._k
    except AttributeError:
        return obj.pyobject()._k


def canonicalize_elem_syms(expr):
    expr = sympify(expr)
    if not expr.args:
        return expr
    expr = expand(expr)
    if is_of_func_type(expr, FactorialElemSym):
        if degree(expr) < numvars(expr):
            return canonicalize_elem_syms(split_out_vars(expr, genvars(expr)[:-1], None))
        return expr
    if isinstance(expr, Add):
        return Add(*[canonicalize_elem_syms(arg) for arg in expr.args])
    if isinstance(expr, Mul):
        # make the z variables disjoint
        if any(isinstance(arg, Add) for arg in expr.args):
            return canonicalize_elem_syms(expand(expr))
        split_out = [arg for arg in expr.args if not is_of_func_type(arg, FactorialElemSym) and not isinstance(arg, Pow)]
        elems = [arg for arg in expr.args if is_of_func_type(arg, FactorialElemSym)]
        pows = [arg for arg in expr.args if isinstance(arg, Pow)]
        if any(isinstance(arg.args[0], Add) for arg in pows):
            return canonicalize_elem_syms(expand(expr))
        # for arg in pows:
        #     elems += [*(int(arg.args[1])*[arg.args[0]])]
        # split out vars if p != k
        for i, elem in enumerate(elems):
            if isinstance(elem, FactorialElemSym) and elem._p < elem._k:
                elems[i] = elem.split_out_vars(genvars(elem)[: len(genvars(elem)) // 2], None)
                return canonicalize_elem_syms(expand(Mul(*elems, *split_out, *pows)))
        # if we got here, all _p == _k
        var_dict = {}
        for elem in elems:
            var_dict[coeffvars(elem)[0]] = [*var_dict.get(coeffvars(elem)[0], []), *genvars(elem)]
        return Mul(*[*split_out, *[FactorialElemSym(len(v), len(v), v, [k]) for k, v in var_dict.items()], *[Pow(canonicalize_elem_syms(arg.args[0]), arg.args[1]) for arg in pows]])
    return expr


def split_out_vars(expr, vars1, vars2):
    expr = sympify(expr)
    if hasattr(expr, "split_out_vars"):
        try:
            return expr.split_out_vars(vars1, vars2)
        except NotEnoughGeneratorsError:
            return expr
    if hasattr(expr, "pyobject") and hasattr(expr.pyobject(), "split_out_vars"):
        try:
            return sympify(expr.pyobject().split_out_vars(vars1, vars2))
        except NotEnoughGeneratorsError:
            return expr
    if not expr.args:
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
    if isinstance(expr, FactorialElemSym):
        return expr

    if not arg:
        for arg in expr.args:
            expr = elem_sym_unify(expr, arg)
        return expr.func(*[elem_sym_unify(arg) for arg in expr.args])
    expr2 = expr
    # if isinstance(arg, FactorialElemSym):
    #     v1 = Wild("v_1")
    #     v2 = Wild("v_2")
    #     rep_pattern = FactorialElemSym(arg._p, arg._k + 1, [*arg.genvars, v1], [*arg.coeffvars, v2])
    #     pattern = arg + FactorialElemSym(1, 1, [v1], [v2]) * FactorialElemSym(arg._p - 1, arg._k, arg.genvars, [*arg.coeffvars, v2])
    #     expr2 = expr.replace(pattern, rep_pattern)
    if arg.args:
        for arg2 in arg.args:
            expr2 = elem_sym_unify(expr2, arg2)
    if expr != expr2:
        return elem_sym_unify(expr2)
    return expr
