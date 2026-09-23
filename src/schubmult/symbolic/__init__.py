"""Symbolic computation facade: fast SymEngine arithmetic with SymPy printing and polynomial domains.

Import ``Add``, ``Mul``, ``Pow``, ``S``, ``Symbol``, ``sympify``, ``expand`` from here rather than
from ``symengine``/``sympy`` directly; the SymPy versions are available under ``sympy_``-prefixed
names (``sympy_Add``, ``sympy_Mul``, ``sympify_sympy``, ``sympy_poly``). Also re-exports the
Schubert-polynomial helpers of `schubmult.symbolic.poly.schub_poly` and SymPy's
``EXRAW``/``CoercionFailed`` used by the ring domains. Generating sets live in
`schubmult.symbolic.poly.variables` (re-exported from `schubmult.symbolic.poly`).
Exports resolve lazily (PEP 562), so importing a submodule such as
`schubmult.symbolic.poly.variables` does not load SymPy/SymEngine until a name that needs them is used.
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
    "CompositeDomain": ("sympy.polys.domains.compositedomain", "CompositeDomain"),
    "DomainElement": ("sympy.polys.domains.domainelement", "DomainElement"),
    "EXRAW": ("sympy.polys.domains.expressionrawdomain", "EXRAW"),
    "Ring": ("sympy.polys.domains.ring", "Ring"),
    "CoercionFailed": ("sympy.polys.polyerrors", "CoercionFailed"),
    "DefaultPrinting": ("sympy.printing.defaults", "DefaultPrinting"),
}
_sympy.update({name: ("sympy", name) for name in ("Expr", "Function", "Poly", "expand_func", "init_printing", "latex", "poly", "pretty", "prod", "simplify", "sstr")})
_functions = ("efficient_subs", "expand", "expand_seq", "is_of_func_type", "symbols", "sympify")

__all__ = [*_symengine, *_sympy, *_functions]  # noqa: PLE0604


def __getattr__(name):
    import importlib

    if name in _symengine:
        val = getattr(importlib.import_module("symengine"), _symengine[name])
    elif name in _sympy:
        modname, attr = _sympy[name]
        val = getattr(importlib.import_module(modname), attr)
    elif name in _functions:
        val = getattr(importlib.import_module("schubmult.symbolic.functions"), name)
    else:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    globals()[name] = val
    return val


def __dir__():
    return sorted(set(globals()) | set(__all__))
