"""Symbolic computation facade: fast SymEngine arithmetic with SymPy printing and polynomial algorithms.

Import ``Add``, ``Mul``, ``Pow``, ``S``, ``Symbol``, ``sympify``, ``expand`` from here rather than
from ``symengine``/``sympy`` directly; the SymPy versions are available under ``sympy_``-prefixed
names (``sympy_Add``, ``sympy_Mul``, ``sympify_sympy``, ``sympy_poly``). The ring-domain protocol
(``EXRAW``, ``CoercionFailed``, ``Ring``, ...) comes from `schubmult.symbolic.domain`. Generating
sets live in `schubmult.symbolic.poly.variables` (re-exported from `schubmult.symbolic.poly`).

Nothing here imports SymPy eagerly: SymPy names are `schubmult.utils._lazy.LazyAttr` proxies that import on first use
(call, attribute access, or subclassing), so SymPy loads only when something actually needs it.
"""

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from symengine import Add, Integer, Mul, Pow, S, Symbol, SympifyError
    from symengine import sympify as sympify_symengine
    from sympy import (
        Add as sympy_Add,
    )
    from sympy import (
        Expr,
        Function,
        Poly,
        Tuple,
        UnevaluatedExpr,
        expand_func,
        init_printing,
        latex,
        poly,
        pretty,
        pretty_print,
        simplify,
        sstr,
    )
    from sympy import (
        Mul as sympy_Mul,
    )
    from sympy import (
        Symbol as sympy_Symbol,
    )
    from sympy import (
        expand as sympy_expand,
    )
    from sympy import (
        poly as sympy_poly,
    )
    from sympy import (
        sympify as sympify_sympy,
    )
    from sympy.parsing.sympy_parser import parse_expr

    from schubmult.symbolic.domain import (
        EXRAW,
        CoercionFailed,
        CompositeDomain,
        DomainElement,
        Ring,
    )
    from schubmult.symbolic.functions import (
        efficient_subs,
        expand_seq,
        is_of_func_type,
        prod,
        symbols,
        sympify,
    )
    from schubmult.utils._printable import LazyPrintable as DefaultPrinting

    def expand(obj: Any, **kwargs: Any) -> Any: ...


_symengine = {
    "Add": "Add",
    "Integer": "Integer",
    "Mul": "Mul",
    "Pow": "Pow",
    "S": "S",
    "Symbol": "Symbol",
    "SympifyError": "SympifyError",
    "sympify_symengine": "sympify",
}
_sympy = {
    "sympy_Add": ("sympy", "Add"),
    "sympy_Mul": ("sympy", "Mul"),
    "sympify_sympy": ("sympy", "sympify"),
    "sympy_poly": ("sympy", "poly"),
    "sympy_expand": ("sympy", "expand"),
    "sympy_Symbol": ("sympy", "Symbol"),
    "parse_expr": ("sympy.parsing.sympy_parser", "parse_expr"),
}
_sympy.update(
    {name: ("sympy", name) for name in ("Expr", "Function", "Poly", "Tuple", "UnevaluatedExpr", "expand_func", "init_printing", "latex", "poly", "pretty", "pretty_print", "simplify", "sstr")},
)
_domain = ("CoercionFailed", "CompositeDomain", "DomainElement", "EXRAW", "Ring")
_functions = ("efficient_subs", "expand", "expand_seq", "is_of_func_type", "prod", "symbols", "sympify")

__all__ = [
    "EXRAW",
    "Add",
    "CoercionFailed",
    "CompositeDomain",
    "DefaultPrinting",
    "DomainElement",
    "Expr",
    "Function",
    "Integer",
    "Mul",
    "Poly",
    "Pow",
    "Ring",
    "S",
    "Symbol",
    "SympifyError",
    "Tuple",
    "UnevaluatedExpr",
    "efficient_subs",
    "expand",
    "expand_func",
    "expand_seq",
    "init_printing",
    "is_of_func_type",
    "latex",
    "parse_expr",
    "poly",
    "pretty",
    "pretty_print",
    "prod",
    "simplify",
    "sstr",
    "symbols",
    "sympify",
    "sympify_symengine",
    "sympify_sympy",
    "sympy_Add",
    "sympy_Mul",
    "sympy_Symbol",
    "sympy_expand",
    "sympy_poly",
]


def __getattr__(name: str) -> Any:
    import importlib

    if name in _symengine:
        val = getattr(importlib.import_module("symengine"), _symengine[name])
    elif name in _sympy:
        from schubmult.utils._lazy import LazyAttr

        val = LazyAttr(*_sympy[name])
    elif name in _domain:
        val = getattr(importlib.import_module("schubmult.symbolic.domain"), name)
    elif name in _functions:
        val = getattr(importlib.import_module("schubmult.symbolic.functions"), name)
    elif name == "DefaultPrinting":
        from schubmult.utils._printable import LazyPrintable as val
    else:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    globals()[name] = val
    return val


def __dir__():
    return sorted(set(globals()) | set(__all__))
