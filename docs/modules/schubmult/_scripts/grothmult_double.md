<a id="schubmult._scripts.grothmult_double"></a>

# schubmult.\_scripts.grothmult\_double

``grothmult_double`` console script: products of double Grothendieck polynomials.

**Example**:

  grothmult_double 3 1 2 - 2 1 3 --display-positive
  grothmult_double --code 2 0 - 1 0 --mixed-var

<a id="schubmult._scripts.grothmult_double.groth_posify"></a>

#### groth\_posify

```python
def groth_posify(val, var2, var3, msg)
```

Positive FGL form of a coefficient, at ``beta = 1``.

Basis: all products of ``z_i - y_j``, ``(1 + y_j)^{-1}``, ``(1 + z_i)^{-1}``.
Substituting ``u_j = 1 + y_j``, ``v_i = 1 + z_i`` (so ``y = u - 1``,
``z = v - 1``) turns the value into a Laurent polynomial in ``u, v``; clearing
the monomial denominator makes everything polynomial.  Candidates are
difference products times leftover clearing monomials, expanded once, with
``as_coefficients_dict`` terms used as opaque basis vectors for an integer
LP.  ``beta`` is restored as ``beta**(`diffs` - d)`` with the Laurent atoms
``(1 + beta*y)`` degree 0 by construction.

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
``--display-positive`` (only supported together with ``--mixed-var``), coefficients
are rewritten in the positive FGL basis of ``z_i - y_j``, ``(1 + y_j)^{-1}``, and
``(1 + z_i)^{-1}`` via an integer program (see `groth_posify`). Returns the raw
result dict when the caller passes a ``None`` formatter (e.g. from tests); otherwise
prints and returns ``None``.

