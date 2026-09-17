<a id="schubmult.mult.groth"></a>

# schubmult.mult.groth

beta-Grothendieck Chevalley formula: multiplication of a (single) Grothendieck
polynomial by a bare x_k variable, i.e. x_k * G_w^(beta).

Non-equivariant (y=0) for now; ``GrothendieckRing`` has no coefficient/y genset yet.

Derived from M. Willems, "A Chevalley formula in equivariant K-theory"
(arXiv:math/0603220), Theorem 5 (the ordinary, non-equivariant specialization of
his equivariant Chevalley formula, Theorem 4). Willems indexes K-theory classes
O_w by the *dimension* of the Schubert variety, dual to the *codimension*
indexing used by Schubert/Grothendieck polynomials S_w/G_w; the w0-conjugation
below (``hat_w = w0*w`` going in, ``w0*v`` coming out) translates between the two
conventions. The beta-grading (beta^(d-1) per length difference d = l(v)-l(w))
matches this codebase's beta-deformed Grothendieck polynomial normalization
(beta=0 recovers the classical double Schubert Monk formula). Calibrated against
grothendieck_poly()/to_groth() (see session notes).

<a id="schubmult.mult.groth.chevalley_x_k"></a>

#### chevalley\_x\_k

```python
def chevalley_x_k(w, k, beta, n=None)
```

Coefficients of ``x_k * G_w^(beta)`` in the Grothendieck basis, as a dict
``{v: coeff}`` (``w`` itself never appears: the self-term cancels identically).

<a id="schubmult.mult.groth.single_variable_groth"></a>

#### single\_variable\_groth

```python
def single_variable_groth(coeff_dict, varnum, beta)
```

Multiply ``sum_u coeff_u G_u^(beta)`` by the single variable ``x_varnum``
(Grothendieck Chevalley formula), via ``chevalley_x_k``.

<a id="schubmult.mult.groth.mult_poly_groth"></a>

#### mult\_poly\_groth

```python
def mult_poly_groth(coeff_dict, poly, var_x, beta)
```

Multiply ``sum_u coeff_u G_u^(beta)`` by an arbitrary polynomial ``poly`` in ``var_x``.

Recurses over the ``Add``/``Mul``/``Pow`` structure of ``poly``, dispatching
single-variable leaves to ``single_variable_groth``.

