"""Symbolic computation facade: fast SymEngine arithmetic with SymPy printing and polynomial domains.

Import ``Add``, ``Mul``, ``Pow``, ``S``, ``Symbol``, ``sympify``, ``expand`` from here rather than
from ``symengine``/``sympy`` directly; the SymPy versions are available under ``sympy_``-prefixed
names (``sympy_Add``, ``sympy_Mul``, ``sympify_sympy``, ``sympy_poly``). Also re-exports the
Schubert-polynomial helpers of `schubmult.symbolic.poly.schub_poly` and SymPy's
``EXRAW``/``CoercionFailed`` used by the ring domains. Generating sets live in
`schubmult.symbolic.poly.variables` (re-exported from `schubmult.symbolic.poly`).
"""

# sympy.init_printing(str_printer = sstr)
from symengine import Add, Integer, Mul, Pow, S, Symbol, SympifyError
from symengine import sympify as sympify_symengine
from sympy import Add as sympy_Add
from sympy import Expr, Function, Poly, expand_func, init_printing, latex, poly, pretty, prod, simplify, sstr
from sympy import Mul as sympy_Mul
from sympy import sympify as sympify_sympy
from sympy.polys.domains.compositedomain import CompositeDomain
from sympy.polys.domains.domainelement import DomainElement
from sympy.polys.domains.expressionrawdomain import EXRAW
from sympy.polys.domains.ring import Ring
from sympy.polys.polyerrors import CoercionFailed
from sympy.printing.defaults import DefaultPrinting

from .functions import efficient_subs, expand, expand_seq, is_of_func_type, symbols, sympify

sympy_poly = poly
