<a id="schubmult._scripts.schubmult_py"></a>

# schubmult.\_scripts.schubmult\_py

``schubmult_py`` console script: products of ordinary Schubert polynomials.

**Example**:

  schubmult_py 3 1 2 - 2 1 3
  schubmult_py --code 2 0 - 1 0
  schubmult_py --coprod --code 2 0 3 0 1 - 2 4

<a id="schubmult._scripts.schubmult_py.main"></a>

#### main

```python
def main(argv=None)
```

Entry point for the ``schubmult_py`` console script.

Parses permutations (or, with ``--code``, Lehmer codes) from ``argv``, multiplies
their ordinary Schubert polynomials via `schubmult.mult.single.schubmult_py`, and
prints the resulting coefficient dictionary ``{Permutation: coefficient}``. With
``--coprod``, computes a coproduct instead (see `schub_coprod_py`). Returns the raw
result dict when the caller passes a ``None`` formatter (e.g. from tests); otherwise
prints and returns ``None``.

