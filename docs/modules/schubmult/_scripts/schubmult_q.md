<a id="schubmult._scripts.schubmult_q"></a>

# schubmult.\_scripts.schubmult\_q

``schubmult_q`` console script: products of quantum Schubert polynomials.

**Example**:

  schubmult_q 3 1 2 - 2 1 3
  schubmult_q 3 1 2 - 2 1 3 --parabolic 1

<a id="schubmult._scripts.schubmult_q.main"></a>

#### main

```python
def main(argv=None)
```

Entry point for the ``schubmult_q`` console script.

Parses permutations (or, with ``--code``, Lehmer codes) from ``argv``, multiplies
their quantum Schubert polynomials via `schubmult.mult.quantum.schubmult_q_fast`
(or the unoptimized `schubmult.mult.quantum.schubmult_q` with ``--basic-pieri``),
and prints the resulting coefficient dictionary ``{Permutation: coefficient}``
(coefficients are polynomials in the quantum parameters ``q_i``). With
``--parabolic g1 g2 ...`` (generator block sizes of a parabolic subgroup, exactly
two input permutations required), the result is projected via the Peterson-Woodward
theorem (`schubmult.mult.quantum_double.apply_peterson_woodward`). Returns the raw
result dict when the caller passes a ``None`` formatter (e.g. from tests); otherwise
prints and returns ``None``.

