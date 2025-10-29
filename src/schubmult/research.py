# metaclass
# symengine and sympy together but otherwise itself
import symengine.lib.symengine_wrapper as sw
import sympy


class Bacon(sympy.Expr, metaclass=sw.BasicMeta):
    pass
