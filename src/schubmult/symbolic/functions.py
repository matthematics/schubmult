import symengine
import symengine.lib.symengine_wrapper as sw
import sympy


def expand(obj, **kwargs):
    if len(kwargs.keys()):
        return symengine.sympify(sympy.expand(obj, **kwargs))
    try:
        return symengine.expand(obj)
    except Exception:
        return sympy.expand(obj)


def symbols(*args, **kwargs):
    return symengine.symbols(*args, **kwargs)


def sympify(val):
    try:
        return symengine.sympify(val)
    except symengine.SympifyError:
        return sympy.sympify(val)


def is_of_func_type(elem, typ):
    return isinstance(elem, typ) or (isinstance(elem, sw.PyFunction) and isinstance(elem.pyobject(), typ))


def expand_seq(seq, genset):
    return sympy.prod([genset[i + 1] ** seq[i] for i in range(len(seq))])
