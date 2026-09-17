<a id="schubmult._scripts.quasi_dd_test"></a>

# schubmult.\_scripts.quasi\_dd\_test

<a id="schubmult._scripts.quasi_dd_test.t_index_iA"></a>

#### t\_index\_iA

```python
def t_index_iA(i, A)
```

Return j such that t_{i, A} = t_j, i.e. j = (\bar A)_i (1-indexed i-th elt of N \ A).

<a id="schubmult._scripts.quasi_dd_test.ev_A"></a>

#### ev\_A

```python
def ev_A(pol, x_arr, A, t_gen, n_xvars)
```

ev_A f = f(t_{(bar A)_1}, ..., t_{(bar A)_{n_xvars}}; t) — substitute x_j -> t_{(bar A)_j}.

<a id="schubmult._scripts.quasi_dd_test.factorization_from_code"></a>

#### factorization\_from\_code

```python
def factorization_from_code(code)
```

Columns factorization of forest F=1^{c_1}·2^{c_2}·... in Thompson monoid.

<a id="schubmult._scripts.quasi_dd_test.a_F_polynomial"></a>

#### a\_F\_polynomial

```python
def a_F_polynomial(f_poly, x_arr, t_gen, code, n_xvars)
```

Compute a_F(t) = [ev ⋆ E_F] f_poly  for forest with given code via direct
polynomial divided differences.

For factorization (i_1, ..., i_k) of F:
    a_F = ev_{A_0} ∘ E_{i_1, A_1} ∘ ... ∘ E_{i_k, ∅}  f
where A_j = i_{j+1} ⋆ ... ⋆ i_k.

<a id="schubmult._scripts.quasi_dd_test.schub_a_F_polynomial"></a>

#### schub\_a\_F\_polynomial

```python
def schub_a_F_polynomial(schub, x_arr, t_gen, code, n_xvars)
```

Compute a_F(t) = [ev ⋆ E_F] f_poly  for forest with given code via direct
polynomial divided differences.

For factorization (i_1, ..., i_k) of F:
    a_F = ev_{A_0} ∘ E_{i_1, A_1} ∘ ... ∘ E_{i_k, ∅}  f
where A_j = i_{j+1} ⋆ ... ⋆ i_k.

<a id="schubmult._scripts.quasi_dd_test.lrcoeff_a_F_polynomial"></a>

#### lrcoeff\_a\_F\_polynomial

```python
def lrcoeff_a_F_polynomial(schub, x_arr, t_gen, code, n_xvars)
```

Compute a_F(t) = [ev ⋆ E_F] f_poly  for forest with given code via direct
polynomial divided differences.

For factorization (i_1, ..., i_k) of F:
    a_F = ev_{A_0} ∘ E_{i_1, A_1} ∘ ... ∘ E_{i_k, ∅}  f
where A_j = i_{j+1} ⋆ ... ⋆ i_k.

<a id="schubmult._scripts.quasi_dd_test.enum_forest_codes"></a>

#### enum\_forest\_codes

```python
def enum_forest_codes(length, max_sum)
```

All forest codes (c_1, ..., c_length) with sum <= max_sum, c_i >= 0.

