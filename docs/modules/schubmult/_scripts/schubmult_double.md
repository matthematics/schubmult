<a id="schubmult._scripts.schubmult_double"></a>

# schubmult.\_scripts.schubmult\_double

``schubmult_double`` console script: products of double Schubert polynomials.

**Example**:

  schubmult_double 3 1 2 - 2 1 3 --display-positive
  schubmult_double --code 2 0 - 1 0 --mixed-var

<a id="schubmult._scripts.schubmult_double.main"></a>

#### main

```python
def main(argv=None)
```

Entry point for the ``schubmult_double`` console script.

Parses permutations (or, with ``--code``, Lehmer codes) from ``argv``, multiplies
their double Schubert polynomials via `schubmult.mult.double.schubmult_double`, and
prints the resulting coefficient dictionary ``{Permutation: coefficient}`` in the
``y``/``z`` coefficient variables. By default the two coefficient variable sets are
the same (``y``); pass ``--mixed-var`` to use distinct ``y``/``z`` sets. With
``--display-positive``, coefficients are rewritten as positive combinations of
``y_i - z_j`` (or verified positive in the single-variable case) via an integer
program (`schubmult.mult.positivity.posify`). Returns the raw result dict when the
caller passes a ``None`` formatter (e.g. from tests); otherwise prints and returns
``None``.

