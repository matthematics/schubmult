"""Parsing of coefficient expressions from the command line."""

import symengine
from latex2sympy2_extended import latex2sympy  # noqa: F401
from sympy.parsing.sympy_parser import parse_expr  # noqa: F401

from schubmult.symbolic.poly.variables import GeneratingSet


def parse_coeff(coeff_str, latex=False):
    """Sympify ``coeff_str`` and map any symbol ``name_i`` to ``GeneratingSet(name)[i]`` so the result
    uses the package's interned variables. LaTeX input is not yet supported (returns ``None``).
    """
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
