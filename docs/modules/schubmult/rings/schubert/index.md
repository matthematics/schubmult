<a id="schubmult.rings.schubert"></a>

# schubmult.rings.schubert

Schubert-family rings: the main user-facing algebra interface.

The most commonly used entry points are the ring *instances*:

- ``Sx``: ordinary (single) Schubert polynomials, ``Sx([3, 1, 2]) * Sx([2, 1, 3])``.
- ``DSx``: double Schubert polynomials (second alphabet ``y``).
- ``Gx`` / ``DGx``: (double) Grothendieck polynomials.
- ``QSx`` / ``QDSx``: quantum (double) Schubert polynomials.
- ``QPSx`` / ``QPDSx``: parabolic quantum (double) Schubert polynomials.

Each instance is an object of the corresponding ``*Ring`` class; calling it with
a permutation (or Lehmer code, or a polynomial expression) yields a ``*Element``
that supports ``+``, ``*``, ``.expand()``, and conversion between bases. All ring
classes derive from `BaseSchubertRing` and dispatch their products to the kernels
in `schubmult.mult`.

Everything here is imported lazily (see ``__getattr__``) so ``import schubmult``
stays fast.

<a id="schubmult.rings.schubert.__getattr__"></a>

#### \_\_getattr\_\_

```python
def __getattr__(name: str)
```

Lazily import and cache the requested export from its defining submodule.

<a id="schubmult.rings.schubert.__dir__"></a>

#### \_\_dir\_\_

```python
def __dir__()
```

Include lazily-exported names in ``dir()``.

