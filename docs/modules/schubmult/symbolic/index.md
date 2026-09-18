<a id="schubmult.symbolic"></a>

# schubmult.symbolic

Symbolic computation facade: fast SymEngine arithmetic with SymPy printing and polynomial domains.

Import ``Add``, ``Mul``, ``Pow``, ``S``, ``Symbol``, ``sympify``, ``expand`` from here rather than
from ``symengine``/``sympy`` directly; the SymPy versions are available under ``sympy_``-prefixed
names (``sympy_Add``, ``sympy_Mul``, ``sympify_sympy``, ``sympy_poly``). Also re-exports the
Schubert-polynomial helpers of `schubmult.symbolic.poly.schub_poly` and SymPy's
``EXRAW``/``CoercionFailed`` used by the ring domains. Generating sets live in
`schubmult.symbolic.poly.variables` (re-exported from `schubmult.symbolic.poly`).

