<a id="schubmult.mult.double"></a>

# schubmult.mult.double

Double Schubert polynomial multiplication.

Implements the ``schubmult_double`` kernel: the product of a linear
combination of double Schubert polynomials ``S_u(x, var2)`` with a single
``S_v(x, var3)``, returned as a coefficient dict ``{w: coeff}`` of polynomials
in ``var2``/``var3``. Uses the same ``theta``/``vmu``/v-path recursion as
``schubmult.mult.single`` (see that module), replacing the plain elementary
symmetric contribution with ``elem_sym_func``, which carries the secondary
variables ``var2``/``var3``.

Also provides the "alt"/"from_elems" variants (building the product one
descent-pulled variable at a time via ``pull_out_var``, generic to any choice
of elementary-symmetric-like function), ``nilhecke_mult`` (nilHecke ring
multiplication), and ``schub_coprod_double`` (the double Schubert coproduct).

<a id="schubmult.mult.double.count_sorted"></a>

#### count\_sorted

```python
def count_sorted(mn, tp)
```

Count occurrences of ``tp`` in the sorted sequence ``mn`` via binary search.

<a id="schubmult.mult.double.single_variable"></a>

#### single\_variable

```python
def single_variable(coeff_dict, varnum, var2=None)
```

Multiply ``sum_u coeff_u S_u(x, var2)`` by the single variable ``x_varnum``.

Equivariant Monk rule: the diagonal term contributes ``var2[u(varnum)]``
(localization of ``x_varnum`` at ``u``) and the off-diagonal terms are the
same Bruhat-cover moves as the ordinary (non-equivariant) ``single_variable``
in ``schubmult.mult.single``.

**Arguments**:

- `coeff_dict` - Mapping ``{Permutation: coeff}``.
- `varnum` - 1-indexed variable index ``k``.
- `var2` - Secondary (``y``) generating set.
  

**Returns**:

- `dict` - The updated coefficient dict.

<a id="schubmult.mult.double.single_variable_down"></a>

#### single\_variable\_down

```python
def single_variable_down(coeff_dict, varnum, var2=None)
```

Down (descent) variant of ``single_variable``, using ``elem_sym_perms_op``.

<a id="schubmult.mult.double.mult_poly_double"></a>

#### mult\_poly\_double

```python
def mult_poly_double(coeff_dict, poly, var_x=None, var_y=None)
```

Multiply ``sum_u coeff_u S_u(x, var_y)`` by an arbitrary polynomial ``poly`` in ``var_x``.

Recurses over the ``Add``/``Mul``/``Pow`` structure of ``poly``, dispatching
single-variable leaves to ``single_variable``; mirrors ``mult_poly_py`` with
the extra secondary alphabet ``var_y``.

<a id="schubmult.mult.double.mult_poly_double_alt"></a>

#### mult\_poly\_double\_alt

```python
def mult_poly_double_alt(coeff_dict, poly, var_x=None, var_y=None)
```

Variant of ``mult_poly_double`` that folds each factor via ``schubmult_double_dict``
instead of ``single_variable``, so ``poly`` is only ever expanded one variable/factor
at a time in the ``S_v`` basis rather than left as a raw scalar multiplier.

<a id="schubmult.mult.double.mult_poly_down"></a>

#### mult\_poly\_down

```python
def mult_poly_down(coeff_dict, poly)
```

Down (descent) variant of ``mult_poly_double``, using ``single_variable_down``
and the fixed default alphabet ``_vars.var1``.

<a id="schubmult.mult.double.nilhecke_mult"></a>

#### nilhecke\_mult

```python
def nilhecke_mult(coeff_dict1, coeff_dict2)
```

NilHecke ring product of ``coeff_dict1`` (polynomial coefficients) and
``coeff_dict2`` (permutation coefficients acting as divided-difference operators).

For each ``w`` in ``coeff_dict2`` its coefficient polynomial is pushed through
``mult_poly_down`` against ``coeff_dict1``, and each resulting permutation ``v``
is right-multiplied by ``w`` whenever that multiplication is length-additive.

**Returns**:

- `dict` - Coefficient dict ``{Permutation: coeff}``.

<a id="schubmult.mult.double.schubmult_double_pair"></a>

#### schubmult\_double\_pair

```python
@cache
def schubmult_double_pair(perm1, perm2, var2=None, var3=None)
```

``schubmult_double`` specialized to a single ``perm1`` with coefficient 1, cached.

<a id="schubmult.mult.double.schubmult_double_pair_generic"></a>

#### schubmult\_double\_pair\_generic

```python
@cache
def schubmult_double_pair_generic(perm1, perm2)
```

``schubmult_double_pair`` with the fixed generic secondary alphabets ``_vars.var_g1``/``_vars.var_g2``.

<a id="schubmult.mult.double.schubmult_double_pair_generic_alt"></a>

#### schubmult\_double\_pair\_generic\_alt

```python
@cache
def schubmult_double_pair_generic_alt(perm1, perm2)
```

Like ``schubmult_double_pair_generic`` but computed via ``schubmult_double_alt_from_elems``
with the factorial elementary symmetric function, then expanded/simplified.

<a id="schubmult.mult.double.schubmult_double_dict"></a>

#### schubmult\_double\_dict

```python
def schubmult_double_dict(perm_dict1, perm_dict2, var2=None, var3=None)
```

Multiply two coefficient dicts of double Schubert polynomials together.

Computes ``(sum_u coeff1_u S_u(x, var2)) * (sum_v coeff2_v S_v(x, var3))``
by summing ``schubmult_double(perm_dict1, v, var2, var3)`` scaled by
``coeff2_v`` over ``v`` in ``perm_dict2``.

<a id="schubmult.mult.double.schubmult_double"></a>

#### schubmult\_double

```python
def schubmult_double(perm_dict, v, var2=None, var3=None)
```

Multiply ``sum_u coeff_u S_u(x, var2)`` by the double Schubert polynomial ``S_v(x, var3)``.

Dispatches to the compiled ``schubmult_cpp`` kernel when available (and both
secondary alphabets are given), falling back to the pure-Python
implementation ``_schubmult_double_python`` otherwise.

**Arguments**:

- `perm_dict` - Mapping ``{Permutation: coeff}``.
- `v` - Permutation (or array-form list) indexing the Schubert polynomial to
  multiply by.
- `var2` - Secondary alphabet attached to ``perm_dict``'s permutations.
- `var3` - Secondary alphabet attached to ``v``.
  

**Returns**:

- `dict` - Coefficient dict ``{Permutation: coeff}`` (polynomials in ``var2``/``var3``).

<a id="schubmult.mult.double.schubmult_double_alt"></a>

#### schubmult\_double\_alt

```python
def schubmult_double_alt(perm_dict, v, var2=None, var3=None, index=1)
```

Alternate double Schubert product, built by peeling one variable of ``~v`` at a
time via ``pull_out_var`` instead of the ``theta``/v-path recursion.

Multiplies ``sum_u coeff_u S_u(x, var2)`` by ``S_v(x, var3)``, recursing on
``~new_v`` with the elementary symmetric factor coming from
``elem_sym_positional_perms`` at each step.

<a id="schubmult.mult.double.schubmult_double_alt_from_elems_forwards"></a>

#### schubmult\_double\_alt\_from\_elems\_forwards

```python
def schubmult_double_alt_from_elems_forwards(perm_dict,
                                             v,
                                             var2=None,
                                             var3=None,
                                             index=1,
                                             elem_func=None)
```

``schubmult_double_alt`` generalized to an arbitrary elementary-symmetric-like
``elem_func(p, k, x_vars, y_vars)``, processing variables of ``~v`` from the first
pulled-out index forward.

<a id="schubmult.mult.double.schubmult_double_alt_from_elems_backwards"></a>

#### schubmult\_double\_alt\_from\_elems\_backwards

```python
def schubmult_double_alt_from_elems_backwards(perm_dict,
                                              v,
                                              var2=None,
                                              var3=None,
                                              elem_func=None)
```

Like ``schubmult_double_alt_from_elems_forwards`` but processing ``~v``'s pulled-out
variables from the last descent backward, multiplying the elementary-symmetric
factor in *before* recursing (dispatches to the compiled kernel when available).

<a id="schubmult.mult.double.schubmult_double_alt_from_elems_backwards_backwards"></a>

#### schubmult\_double\_alt\_from\_elems\_backwards\_backwards

```python
def schubmult_double_alt_from_elems_backwards_backwards(
        perm_dict, v, var2=None, var3=None, elem_func=None)
```

Variant of ``_schubmult_double_alt_from_elems_backwards_python`` without the
per-``new_v`` memoization cache, recursing on the interim dict instead of the
original ``perm_dict`` at each pulled-out variable.

<a id="schubmult.mult.double.schubmult_double_from_elems"></a>

#### schubmult\_double\_from\_elems

```python
def schubmult_double_from_elems(perm_dict,
                                v,
                                var2=None,
                                var3=None,
                                elem_func=None)
```

``schubmult_double`` generalized to an arbitrary elementary-symmetric-like
``elem_func``, via the ``theta``/v-path recursion (rather than ``pull_out_var``).

Dispatches to the compiled kernel when available, falling back to
``_schubmult_double_from_elems_python``.

<a id="schubmult.mult.double.schubmult_double_down"></a>

#### schubmult\_double\_down

```python
def schubmult_double_down(perm_dict, v, var2=None, var3=None)
```

Down (descent) variant of ``_schubmult_double_python``, using ``elem_sym_perms_op``.

<a id="schubmult.mult.double.schub_coprod_double"></a>

#### schub\_coprod\_double

```python
def schub_coprod_double(mperm, indices, var2=None, var3=None)
```

Coproduct of the double Schubert polynomial ``S_mperm`` restricted to the
variable split named by ``indices``.

Analogue of ``schub_coprod_py``: multiplies the Grassmannian permutation for
``indices`` against ``mperm`` (via ``schubmult_double`` with a merged ``2N``
variable alphabet), splits each resulting permutation's window, and
substitutes the merged alphabet back to ``var2``/``var3``.

**Arguments**:

- `mperm` - Permutation (or array-form list) to take the coproduct of.
- `indices` - Iterable of 1-indexed positions selecting the variable split.
- `var2` - Secondary alphabet for the first factor's variables.
- `var3` - Secondary alphabet for the second factor's variables.
  

**Returns**:

- `dict` - Mapping ``{(firstperm, secondperm): coeff}``.

