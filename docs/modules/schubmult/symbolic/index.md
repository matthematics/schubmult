<a id="schubmult.symbolic"></a>

# schubmult.symbolic

Symbolic computation facade: fast SymEngine arithmetic with SymPy printing and polynomial algorithms.

Import ``Add``, ``Mul``, ``Pow``, ``S``, ``Symbol``, ``sympify``, ``expand`` from here rather than
from ``symengine``/``sympy`` directly; the SymPy versions are available under ``sympy_``-prefixed
names (``sympy_Add``, ``sympy_Mul``, ``sympify_sympy``, ``sympy_poly``). The ring-domain protocol
(``EXRAW``, ``CoercionFailed``, ``Ring``, ...) comes from `schubmult.symbolic.domain`. Generating
sets live in `schubmult.symbolic.poly.variables` (re-exported from `schubmult.symbolic.poly`).

Nothing here imports SymPy eagerly: SymPy names are `schubmult.utils._lazy.LazyAttr` proxies that import on first use
(call, attribute access, or subclassing), so SymPy loads only when something actually needs it.

