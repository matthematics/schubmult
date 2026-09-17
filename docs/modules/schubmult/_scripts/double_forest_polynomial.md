<a id="schubmult._scripts.double_forest_polynomial"></a>

# schubmult.\_scripts.double\_forest\_polynomial

Compute double forest polynomials P_F(x; t) via the vine subword model.

Reference
---------
N. Bergeron, L. Gagnon, P. Nadeau, H. Spink, V. Tewari,
"Equivariant quasisymmetry and noncrossing partitions" (arXiv:2504.15234),
Theorem 5.1 and Section 5.1.

Given a code c = (c_1, c_2, ...) for an indexed forest F (so that c(F) = c),
and two sequences of symbols x = (x_1, x_2, ...) and t = (t_1, t_2, ...),
this returns the double forest polynomial

    P_F(x; t) = sum_{pi in R(omega_tilde_[n]; F)} wt(pi),

where omega_tilde_[n] is the long vine word and the sum is over its
subwords whose value-sequence is a Sylvester word for F.

<a id="schubmult._scripts.double_forest_polynomial.forest_from_code"></a>

#### forest\_from\_code

```python
def forest_from_code(code)
```

Build the indexed forest F with c(F) = ``code`` via the Thompson
monoid factorization F = 1^{c_1} . 2^{c_2} . 3^{c_3} . ...

<a id="schubmult._scripts.double_forest_polynomial.sylvester_word"></a>

#### sylvester\_word

```python
def sylvester_word(forest)
```

Return one Sylvester word of ``forest`` via pre-order traversal
(root, left subtree, right subtree) over each non-trivial tree.

<a id="schubmult._scripts.double_forest_polynomial.long_word"></a>

#### long\_word

```python
def long_word(n)
```

Build omega_tilde_[n] as a list of letter records.

omega^(k)_[n] = (n, n-1, ..., k+1, k, k+1_bar, ..., n-1_bar, n_bar)
Concatenate for k = 1, ..., n.

<a id="schubmult._scripts.double_forest_polynomial.letter_weight"></a>

#### letter\_weight

```python
def letter_weight(letter, x_gen, t_gen)
```

wt(j^(i))      = x_i - t_j   (unbarred)
wt(j_bar^(i)) = t_j - t_i   (barred)

<a id="schubmult._scripts.double_forest_polynomial.double_forest_polynomial"></a>

#### double\_forest\_polynomial

```python
def double_forest_polynomial(code, x_gen, t_gen, n=None)
```

Compute P_F(x; t) for the indexed forest F with c(F) = ``code``.

Parameters
----------
code : sequence[int]
    The code (c_1, c_2, ...) of F. Trailing zeros are ignored.
x_gen : callable[int -> Expr]
    Maps i -> x_i  (the non-equivariant variables).
t_gen : callable[int -> Expr]
    Maps j -> t_j  (the equivariant variables).
n : int, optional
    Use the long word omega_tilde_[n]. If ``None``, the smallest n with
    F supported in {1, ..., n+1} is used.

Returns
-------
sympy.Expr
    The expanded polynomial P_F(x; t).

<a id="schubmult._scripts.double_forest_polynomial.decompose_double_forest_tensor"></a>

#### decompose\_double\_forest\_tensor

```python
def decompose_double_forest_tensor(poly, x_genset, t_genset, length)
```

Decompose a double forest polynomial into t-forest ⊗ x-forest terms.

Returns a dict mapping ``(t_code, x_code) -> integer/symbolic coefficient``.

<a id="schubmult._scripts.double_forest_polynomial.DoubleForestPolynomialBasis"></a>

## DoubleForestPolynomialBasis Objects

```python
class DoubleForestPolynomialBasis()
```

Concrete basis for double forest polynomials.

A basis element is indexed by a pair ``(t_code, x_code)`` and represents
``Forest_t(t_code) ⊗ Forest_x(x_code)``.

<a id="schubmult._scripts.double_forest_polynomial.ForestDoubleElement"></a>

## ForestDoubleElement Objects

```python
class ForestDoubleElement(dict)
```

Formal element in the abstract double-forest basis.

Stored as {forest_code: coeff}, representing
    sum_F coeff[F] * DF(F),
where DF(F) is an abstract basis symbol (not expanded into x/t variables).

<a id="schubmult._scripts.double_forest_polynomial.ForestDoubleElement.expand_polynomial"></a>

#### expand\_polynomial

```python
def expand_polynomial()
```

Expand abstract basis expression to a concrete x/t polynomial.

<a id="schubmult._scripts.double_forest_polynomial.ForestDouble"></a>

## ForestDouble Objects

```python
class ForestDouble()
```

Abstract double-forest basis algebra.

Basis symbols are indexed by one forest code F and represent P_F(x; t)
abstractly. Coefficients live in the polynomial ring of equivariant vars.

<a id="schubmult._scripts.double_forest_polynomial.extract_double_forest_coefficients"></a>

#### extract\_double\_forest\_coefficients

```python
def extract_double_forest_coefficients(poly, x_genset, length)
```

Return all coefficients a_F(t) in f = sum_F a_F(t) P_F(x; t).

Implemented by iterative top-degree peeling:
1) convert the current residual into ForestPolyBasis in x,
2) read top-degree coefficients,
3) subtract those coefficients times the corresponding double forest basis
   polynomials P_F,
4) repeat until the residual vanishes.

This avoids treating the double-forest basis as a plain one-shot x-basis
conversion and follows the triangular subtraction strategy.

<a id="schubmult._scripts.double_forest_polynomial.extract_double_forest_coefficient"></a>

#### extract\_double\_forest\_coefficient

```python
def extract_double_forest_coefficient(poly, forest_code, x_genset, length)
```

Return the single coefficient a_F(t) of P_F in f.

Equivalent to the operator formula a_F(t) = [ev star e_F] f from
Theorem 10.9 (Coefficient extraction and star-composition).

Note: this returns a t-polynomial coefficient in the x-forest basis.

<a id="schubmult._scripts.double_forest_polynomial.extract_double_forest_tensor_coefficients"></a>

#### extract\_double\_forest\_tensor\_coefficients

```python
def extract_double_forest_tensor_coefficients(poly, x_genset, t_genset,
                                              length)
```

Return true double-forest tensor-basis coefficients.

Expands poly as
    poly = sum_{A,B} c_{A,B} Forest_t(A) \otimes Forest_x(B),
returning a dict mapping (A, B) -> c_{A,B}.

<a id="schubmult._scripts.double_forest_polynomial.extract_double_forest_tensor_coefficient"></a>

#### extract\_double\_forest\_tensor\_coefficient

```python
def extract_double_forest_tensor_coefficient(poly, t_code, x_code, x_genset,
                                             t_genset, length)
```

Return c_{A,B} for one tensor basis pair A=t_code, B=x_code.

<a id="schubmult._scripts.double_forest_polynomial.monk_style_degree_one_product"></a>

#### monk\_style\_degree\_one\_product

```python
def monk_style_degree_one_product(code,
                                  i,
                                  x_gen,
                                  t_gen,
                                  x_genset,
                                  t_genset,
                                  n=None,
                                  length=None)
```

Compute a Monk-style decomposition for degree-1 double forest multiplication.

The product
    P_{e_i}(x; t) * P_code(x; t)
is decomposed as
    diagonal_coeff(t) * P_code(x; t) + sum_{beta != code} c_beta(t) * P_beta(x; t),
and also expanded in the tensor basis t-forest ⊗ x-forest.

