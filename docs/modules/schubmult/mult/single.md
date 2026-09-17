<a id="schubmult.mult.single"></a>

# schubmult.mult.single

Ordinary (single) Schubert polynomial multiplication.

Implements the ``schubmult_py`` kernel: the product of a linear combination of
Schubert polynomials ``S_u`` (given as a dict ``{u: coeff}``) with a single
Schubert polynomial ``S_v``, returned as a coefficient dict ``{w: coeff}``.

The algorithm is the recursive "v-path" / transition method used throughout
``schubmult``: write ``theta = (~v).theta()`` for the dominant weakly-decreasing
vector bounding ``v``'s Lehmer code, let ``mu = uncode(theta)`` and
``vmu = v * mu``, and process the entries of ``theta`` one at a time, tracking
Bruhat-chain "v-paths" from ``vmu`` down to the identity (``compute_vpathdicts``)
alongside chains of elementary-symmetric moves on the ``u`` side
(``elem_sym_perms``). Accumulating consistent pairs of chains and reading off
the coefficient landing on ``vmu`` gives the product.

<a id="schubmult.mult.single.single_variable"></a>

#### single\_variable

```python
def single_variable(coeff_dict, varnum)
```

Multiply ``sum_u coeff_u S_u(x)`` by the single variable ``x_varnum``.

Uses the classical Monk rule: ``x_k * S_u = sum S_{u t_{ij}}`` over Bruhat
covers ``u t_{ij}`` with ``i <= k < j`` (added) minus those with ``j <= k < i``
(subtracted), via ``elem_sym_perms(u, 1, varnum)``.

**Arguments**:

- `coeff_dict` - Mapping ``{Permutation: coeff}``.
- `varnum` - 1-indexed variable index ``k``.
  

**Returns**:

- `dict` - The updated coefficient dict ``{Permutation: coeff}``.

<a id="schubmult.mult.single.mult_poly_py"></a>

#### mult\_poly\_py

```python
def mult_poly_py(coeff_dict, poly, var_x=_vars.var_x)
```

Multiply ``sum_u coeff_u S_u(x)`` by an arbitrary polynomial ``poly`` in ``var_x``.

Recurses over the ``Add``/``Mul``/``Pow`` structure of ``poly``; each single
variable leaf is dispatched to ``single_variable``, and any other leaf just
scales every coefficient.

**Arguments**:

- `coeff_dict` - Mapping ``{Permutation: coeff}``.
- `poly` - Symbolic polynomial expression in the variables of ``var_x``.
- `var_x` - Generating set identifying the ``x`` variables (default ``x``).
  

**Returns**:

- `dict` - The updated coefficient dict ``{Permutation: coeff}``.

<a id="schubmult.mult.single.schubmult_py"></a>

#### schubmult\_py

```python
def schubmult_py(perm_dict, v)
```

Multiply ``sum_u coeff_u S_u(x)`` by the (ordinary) Schubert polynomial ``S_v``.

Dispatches to the compiled ``schubmult_cpp`` kernel when available and the
permutations fit within its ``MAXN``, falling back to the pure-Python
implementation otherwise.

**Arguments**:

- `perm_dict` - Mapping ``{Permutation: coeff}`` (integer coefficients).
- `v` - Permutation (or array-form list) indexing the Schubert polynomial to
  multiply by.
  

**Returns**:

- `dict` - Coefficient dict ``{Permutation: coeff}`` for the product.

<a id="schubmult.mult.single.schubmult_py_down"></a>

#### schubmult\_py\_down

```python
def schubmult_py_down(perm_dict, v)
```

Divided-difference ("down") variant of ``_schubmult_py_python``.

Same v-path recursion but built from ``elem_sym_perms_op`` (Bruhat *descents*)
instead of ``elem_sym_perms``, used for the down/dual side of the transition
recursion rather than ordinary multiplication.

**Arguments**:

- `perm_dict` - Mapping ``{Permutation: coeff}``.
- `v` - Permutation to multiply by.
  

**Returns**:

- `dict` - Coefficient dict ``{Permutation: coeff}``.

<a id="schubmult.mult.single.schub_coprod_py"></a>

#### schub\_coprod\_py

```python
def schub_coprod_py(perm, indices)
```

Coproduct of ``S_perm`` restricted to the variable split named by ``indices``.

Computes the expansion of the (single) Schubert polynomial coproduct
``Delta_{indices}(S_perm) = sum (firstperm, secondperm) -> coeff`` by
multiplying the Grassmannian permutation for ``indices`` against ``perm``
via ``schubmult_py`` and splitting each resulting permutation's window into
its first ``N`` and remaining ``len(perm) - N`` values.

**Arguments**:

- `perm` - Permutation (or array-form list) to take the coproduct of.
- `indices` - Iterable of 1-indexed positions selecting the variable split.
  

**Returns**:

- `dict` - Mapping ``{(firstperm, secondperm): coeff}``.

