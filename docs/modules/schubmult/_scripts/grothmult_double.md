<a id="schubmult._scripts.grothmult_double"></a>

# schubmult.\_scripts.grothmult\_double

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

