<a id="schubmult._scripts.grothmult_double"></a>

# schubmult.\_scripts.grothmult\_double

``grothmult_double`` console script: products of double Grothendieck polynomials.

**Example**:

  grothmult_double 3 1 2 - 2 1 3 --display-positive
  grothmult_double --code 2 0 - 1 0 --mixed-var

<a id="schubmult._scripts.grothmult_double.main"></a>

#### main

```python
def main(argv=None)
```

Entry point for the ``grothmult_double`` console script.

Parses permutations (or, with ``--code``, Lehmer codes) from ``argv``, multiplies
their double Grothendieck polynomials via `schubmult.mult.groth_double.grothmult_double`,
and prints the resulting coefficient dictionary ``{Permutation: coefficient}`` in the
``y``/``z`` coefficient variables and the K-theory parameter ``beta``. With
``--multiplicative-variables``, use ``beta=-1`` and ``1-1/y_i`` coefficient
coordinates (also ``1-1/z_i`` for mixed variables). With
``--display-positive`` (only supported together with ``--mixed-var``), coefficients
are rewritten in the positive FGL basis of ``z_i - y_j``, ``(1 + y_j)^{-1}``, and
``(1 + z_i)^{-1}`` via an integer program (see `schubmult.mult.positivity.groth_posify`);
the integer programs of different coefficients are solved concurrently, up to ``--jobs`` at
a time (default 1; see `schubmult.mult.positivity.groth_posify_many`). Returns the raw
result dict when the caller passes a ``None`` formatter (e.g. from tests); otherwise
prints and returns ``None``.
