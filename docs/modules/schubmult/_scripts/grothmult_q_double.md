<a id="schubmult._scripts.grothmult_q_double"></a>

# schubmult.\_scripts.grothmult\_q\_double

``grothmult_q_double`` console script: products of quantum double Grothendieck polynomials.

Conjectural quantum K-theory Pieri rule; ``--display-positive``, ``--parabolic``, and
``--nil-hecke``/``--nil-hecke-apply`` are not yet supported.

**Example**:

  grothmult_q_double 3 1 2 - 2 1 3 --mixed-var

``--multiplicative-variables`` uses ``beta=-1`` and ``ChangeOfVariables`` with
coefficient coordinates ``1-1/y_i`` (also ``1-1/z_i`` with ``--mixed-var``).
Use ``--simplify`` to cancel the resulting rational coefficients.
Other beta values and ``--display-positive`` conflict with this mode.

<a id="schubmult._scripts.grothmult_q_double.main"></a>

#### main

```python
def main(argv=None)
```

Entry point for the ``grothmult_q_double`` console script.

Parses permutations (or, with ``--code``, Lehmer codes) from ``argv``, multiplies
their quantum double Grothendieck polynomials via
`schubmult.mult.groth_quantum_double.grothmult_q_double`, and prints the resulting
coefficient dictionary ``{Permutation: coefficient}`` in the ``y``/``z`` coefficient
variables and quantum parameters ``q_i``. ``--display-positive``, ``--parabolic``,
``--nil-hecke``, and ``--mult`` are not yet supported and cause an early exit.
Returns the raw result dict when the caller passes a ``None`` formatter (e.g. from
tests); otherwise prints and returns ``None``.
