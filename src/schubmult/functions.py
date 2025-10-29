from schubmult.rings.variables import NotEnoughGeneratorsError
from schubmult.symbolic import Add, Mul, Pow, S, expand, is_of_func_type, prod, sympify

from .elem_sym import FactorialElemSym


def genvars(obj):
    try:
        return obj.genvars
    except AttributeError:
        try:
            return obj.pyobject().genvars
        except AttributeError:
            raise AttributeError("Object or pyobject does have genvars attribute")


def coeffvars(obj):
    try:
        return obj.coeffvars
    except AttributeError:
        try:
            return obj.pyobject().coeffvars
        except AttributeError:
            raise AttributeError("Object or pyobject does have coeffvars attribute")


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


def canonicalize_elem_syms(expr, combine_equal=False):
    expr = sympify(expr)
    expr = expand(expr)
    if not expr.args:
        return expr
    if is_of_func_type(expr, FactorialElemSym):
        if degree(expr) < numvars(expr):
            return canonicalize_elem_syms(split_out_vars(expr, genvars(expr)[: len(genvars(expr)) // 2], None))
        return expr
    if isinstance(expr, Add):
        blaff = Add(*[canonicalize_elem_syms(arg) for arg in expr.args])
    if isinstance(expr, Pow):
        blaff = Pow(canonicalize_elem_syms(expr.args[0]), expr.args[1])
    if isinstance(expr, Mul):
        blaff = Mul(*[canonicalize_elem_syms(arg) for arg in expr.args])
        if blaff == expr:
            mdict = {}
            for arg in blaff.args:
                if not is_of_func_type(arg, FactorialElemSym):
                    mdict[S.One] = mdict.get(S.One, S.One) * arg
                else:
                    check_arg = arg
                    if isinstance(arg, Pow):
                        check_arg = arg.args[0]
                    cv = coeffvars(check_arg)[0]
                    if cv not in mdict:
                        mdict[cv] = arg
                    else:
                        em = mdict[cv]
                        if combine_equal:
                            mdict[cv] = FactorialElemSym(degree(em) + degree(arg), numvars(em) + numvars(arg), *genvars(em), *genvars(arg), cv)
                        else:
                            new_genvars_list = [set()]
                            total_genvars = []
                            if isinstance(mdict[cv], Mul):
                                for arg2 in mdict[cv].args:
                                    total_genvars += genvars(arg2)
                            elif isinstance(mdict[cv], Pow):
                                total_genvars += genvars(mdict[cv].args[0]) * int(mdict[cv].args[1])
                            else:
                                total_genvars += genvars(mdict[cv])
                            total_genvars += genvars(arg)
                            for gv in total_genvars:
                                found = False
                                for L in new_genvars_list:
                                    if gv not in L:
                                        L.add(gv)
                                        found = True
                                        break
                                if not found:
                                    new_genvars_list += [{gv}]
                            # print(new_genvars_list)
                            mdict[cv] = sympify(prod([FactorialElemSym(len(gvs), len(gvs), *gvs, cv) for gvs in new_genvars_list]))
            return Mul(*list(mdict.values()))
    if blaff != expr:
        return canonicalize_elem_syms(blaff)
    return expr


def canonicalize_elem_syms_coeff(expr, combine_equal=False):
    expr = sympify(expr)
    expr = expand(expr)
    # print(f"farfel {expr=}")
    if not expr.args:
        # print(f"I have returned the moofer of sin {expr=}")
        return expr
    if is_of_func_type(expr, FactorialElemSym):
        if degree(expr) < numvars(expr):
            return canonicalize_elem_syms_coeff(split_out_vars(expr, None, coeffvars(expr)[:-1]))
        # print(f"What a bargain {expr=}")
        return expr
    if isinstance(expr, Add):
        blaff = Add(*[canonicalize_elem_syms_coeff(arg) for arg in expr.args])
    if isinstance(expr, Pow):
        blaff = Pow(canonicalize_elem_syms_coeff(expr.args[0]), expr.args[1])
    if isinstance(expr, Mul):
        blaff = Mul(*[canonicalize_elem_syms_coeff(arg) for arg in expr.args])
        if blaff == expr:
            mdict = {}
            for arg in blaff.args:
                if not is_of_func_type(arg, FactorialElemSym):
                    mdict[S.One] = mdict.get(S.One, S.One) * arg
                else:
                    cv = coeffvars(arg)[0]
                    if cv not in mdict:
                        mdict[cv] = arg
                    else:
                        em = mdict[cv]
                        if combine_equal:
                            mdict[cv] = FactorialElemSym(degree(em) + degree(arg), numvars(em) + numvars(arg), *genvars(em), *genvars(arg), cv)
                        else:
                            new_genvars_list = [set()]
                            total_genvars = []
                            if isinstance(mdict[cv], Mul):
                                for arg2 in mdict[cv].args:
                                    total_genvars += genvars(arg2)
                            else:
                                total_genvars += genvars(mdict[cv])
                            total_genvars += genvars(arg)
                            for gv in total_genvars:
                                found = False
                                for L in new_genvars_list:
                                    if gv not in L:
                                        L.add(gv)
                                        found = True
                                        break
                                if not found:
                                    new_genvars_list += [{gv}]
                            # print(new_genvars_list)
                            mdict[cv] = sympify(prod([FactorialElemSym(len(gvs), len(gvs), *gvs, cv) for gvs in new_genvars_list]))
            return Mul(*list(mdict.values()))
    if blaff != expr:
        return canonicalize_elem_syms_coeff(blaff)
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
