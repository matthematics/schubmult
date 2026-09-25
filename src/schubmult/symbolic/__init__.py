"""Symbolic computation facade: fast SymEngine arithmetic with SymPy printing and polynomial algorithms.

Import ``Add``, ``Mul``, ``Pow``, ``S``, ``Symbol``, ``sympify``, ``expand`` from here rather than
from ``symengine``/``sympy`` directly; the SymPy versions are available under ``sympy_``-prefixed
names (``sympy_Add``, ``sympy_Mul``, ``sympify_sympy``, ``sympy_poly``). The ring-domain protocol
(``EXRAW``, ``CoercionFailed``, ``Ring``, ...) comes from `schubmult.symbolic.domain`. Generating
sets live in `schubmult.symbolic.poly.variables` (re-exported from `schubmult.symbolic.poly`).

Nothing here imports SymPy eagerly: SymPy names are `schubmult.utils._lazy.LazyAttr` proxies that import on first use
(call, attribute access, or subclassing), so SymPy loads only when something actually needs it.
"""


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

__all__ = [*_symengine, *_sympy, *_domain, *_functions, "DefaultPrinting"]  # noqa: PLE0604


def __getattr__(name):
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
