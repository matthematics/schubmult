<a id="schubmult._scripts.schubmult_q_double"></a>

# schubmult.\_scripts.schubmult\_q\_double

``schubmult_q_double`` console script: products of quantum double Schubert polynomials.

Conjectural for most cases; see `schubmult.mult.quantum_double`.

**Example**:

  schubmult_q_double 3 1 2 - 2 1 3 --display-positive
  schubmult_q_double 3 1 2 - 2 1 3 --parabolic 1

<a id="schubmult._scripts.schubmult_q_double.main"></a>

#### main

```python
def main(argv=None)
```

Entry point for the ``schubmult_q_double`` console script.

Parses permutations (or, with ``--code``, Lehmer codes) from ``argv``, multiplies
their quantum double Schubert polynomials via
`schubmult.mult.quantum_double.schubmult_q_double_fast` (or the unoptimized
`schubmult.mult.quantum_double.schubmult_q_double` with ``--basic-pieri``), and
prints the resulting coefficient dictionary ``{Permutation: coefficient}`` in the
``y``/``z`` coefficient variables and quantum parameters ``q_i``. With
``--parabolic``, the result is projected via the Peterson-Woodward theorem
(`apply_peterson_woodward`). With ``--display-positive``, coefficients are rewritten
positively via `schubmult.mult.quantum_double.q_posify`. ``--nil-hecke N`` substitutes
up to ``N`` Fomin-Gelfand-Postnikov commuting difference operators
(`schubmult.mult.quantum_double.nil_hecke`) instead of multiplying permutations.
Returns the raw result dict when the caller passes a ``None`` formatter (e.g. from
tests); otherwise prints and returns ``None``.

