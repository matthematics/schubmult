<a id="schubmult._scripts.grothmult_q"></a>

# schubmult.\_scripts.grothmult\_q

``grothmult_q`` console script: products of quantum Grothendieck polynomials.

Conjectural quantum K-theory Pieri rule; ``--parabolic`` and ``--mult`` are not yet
supported.

**Example**:

  grothmult_q 3 1 2 - 2 1 3

<a id="schubmult._scripts.grothmult_q.main"></a>

#### main

```python
def main(argv=None)
```

Entry point for the ``grothmult_q`` console script.

Parses permutations (or, with ``--code``, Lehmer codes) from ``argv``, multiplies
their quantum Grothendieck polynomials via `schubmult.mult.groth_quantum.grothmult_q`,
and prints the resulting coefficient dictionary ``{Permutation: coefficient}``
(coefficients are polynomials in the quantum parameters ``q_i`` and the K-theory
parameter ``beta``). ``--parabolic`` and ``--mult`` are not yet supported and cause
an early exit. Returns the raw result dict when the caller passes a ``None``
formatter (e.g. from tests); otherwise prints and returns ``None``.

