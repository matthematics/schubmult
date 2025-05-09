import symengine
from latex2sympy2_extended import latex2sympy  # noqa: F401
from sympy.parsing.sympy_parser import parse_expr  # noqa: F401

from schubmult.rings.variables import GeneratingSet


def parse_coeff(coeff_str, latex=False):
    if not latex:
        result = symengine.sympify(coeff_str)
        subs_dict = {}
        for s in result.free_symbols:
            if s.name.find("_") != -1:
                base, index = s.name.split("_")
                gset = GeneratingSet(base)
                subs_dict[s] = gset[int(index)]
        return result.subs(subs_dict)
    return None
