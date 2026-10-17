<a id="schubmult._scripts.grothmult_py"></a>

# schubmult.\_scripts.grothmult\_py

``grothmult_py`` console script: products of (ordinary) Grothendieck polynomials.

**Example**:

  grothmult_py 3 1 2 - 2 1 3
  grothmult_py --code 2 0 - 1 0

<a id="schubmult._scripts.grothmult_py.main"></a>

#### main

```python
def main(argv=None)
```

Entry point for the ``grothmult_py`` console script.

Parses permutations (or, with ``--code``, Lehmer codes) from ``argv``, multiplies
their Grothendieck polynomials via `schubmult.mult.groth.grothmult_py`, and prints
the resulting coefficient dictionary ``{Permutation: coefficient}`` (coefficients are
Laurent polynomials in the K-theory parameter ``beta``). Returns the raw result dict
when the caller passes a ``None`` formatter (e.g. from tests); otherwise prints and
returns ``None``.

