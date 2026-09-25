"""Ready-made symbols and generating sets, in the spirit of ``sympy.abc``.

Importing this module creates the SymEngine symbols ``x_1..x_99``, ``y_1..y_99``, ``z_1..z_99``,
``q_1..q_99`` and ``β`` and exposes the following names:

```python
from schubmult.abc import x, y, z, q, beta
```

- ``x``, ``y``, ``z``, ``q`` -- `GeneratingSet` objects; ``x[i]`` is the symbol ``x_i``. ``x`` is
  the default Schubert variable set, ``y``/``z`` the coefficient (double) variable sets, ``q``
  the quantum parameters.
- ``beta`` -- the symbol ``β``, the K-theoretic deformation parameter of Grothendieck polynomials.

It also re-exports the unevaluated symmetric polynomial atoms of
`schubmult.symbolic.symmetric_polynomials`:

```python
from schubmult.abc import E, e, H, h, E_q
```

- ``E`` (alias `FactorialElemSym`) and ``e`` (alias `ElemSym`) -- factorial and ordinary
  elementary symmetric polynomials ``E(p, k, xvars, yvars)`` / ``e(p, k, xvars)``.
- ``H`` (alias `FactorialCompleteSym`) and ``h`` (alias `CompleteSym`) -- the complete
  homogeneous analogues.
- ``E_q`` (alias `QFactorialElemSym`) -- quantum factorial elementary symmetric polynomials.

Example:

```python
>>> from schubmult import Sx
>>> from schubmult.abc import x, y, E
>>> Sx([3, 1, 2]).expand()
x_1**2
>>> E(1, 2, [x[1], x[2]], [y[1], y[2]]).expand(func=True)
x_1 + x_2 - y_1 - y_2
```
"""

from symengine import var

from schubmult.symbolic.poly.variables import GeneratingSet

_symmetric_names = ("CompleteSym", "E", "E_q", "ElemSym", "FactorialCompleteSym", "FactorialElemSym", "H", "QFactorialElemSym", "e", "h")


def __getattr__(name):
    # The symmetric polynomial atoms are SymPy Functions; import them only when requested
    if name in _symmetric_names:
        import schubmult.symbolic.symmetric_polynomials as symmetric_polynomials

        val = getattr(symmetric_polynomials, name)
        globals()[name] = val
        return val
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


var("x_(1:100)")
var("y_(1:100)")
var("z_(1:100)")
var("q_(1:100)")
var("β")

x = GeneratingSet("x")
y = GeneratingSet("y")
z = GeneratingSet("z")
q = GeneratingSet("q")

beta = var("β")
