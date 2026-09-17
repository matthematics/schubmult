<a id="schubmult.mult.quantum"></a>

# schubmult.mult.quantum

Quantum (single) Schubert polynomial multiplication.

Implements ``schubmult_q``/``schubmult_q_fast``: the product of a linear
combination of quantum Schubert polynomials ``S_u`` with a single ``S_v``,
returned as a coefficient dict ``{w: coeff}`` polynomial in the quantum
parameters ``q_1, q_2, ...``. Uses the same ``theta``/v-path recursion as
``schubmult.mult.single``, with the elementary-symmetric step generalized to
the quantum Pieri-type moves of ``elem_sym_perms_q`` (which may pick up a
factor of ``q`` when a Bruhat move is replaced by its quantum analogue).
``schubmult_q`` uses ``strict_theta`` (no repeated layer merging);
``schubmult_q_fast``/``_schubmult_q_fast_python`` uses ``medium_theta`` and
merges adjacent equal-length layers via ``double_elem_sym_q`` for speed.

<a id="schubmult.mult.quantum.single_variable"></a>

#### single\_variable

```python
def single_variable(coeff_dict, varnum, var_q=_vars.q_var)
```

Multiply ``sum_u coeff_u S_u(x)`` by the single variable ``x_varnum`` (quantum Monk rule).

Same structure as ``schubmult.mult.single.single_variable``, using
``elem_sym_perms_q`` so that some Bruhat moves carry a factor from ``var_q``.

<a id="schubmult.mult.quantum.mult_poly_q"></a>

#### mult\_poly\_q

```python
def mult_poly_q(coeff_dict, poly, var_x=_vars.var_x, var_q=_vars.q_var)
```

Multiply ``sum_u coeff_u S_u(x)`` by an arbitrary polynomial ``poly`` in ``var_x``.

Recurses over the ``Add``/``Mul``/``Pow`` structure of ``poly``, dispatching
single-variable leaves to ``single_variable``; mirrors ``mult_poly_py``.

<a id="schubmult.mult.quantum.schubmult_q_fast"></a>

#### schubmult\_q\_fast

```python
def schubmult_q_fast(perm_dict, v, q_var=_vars.q_var)
```

Multiply ``sum_u coeff_u S_u(x)`` by the quantum Schubert polynomial ``S_v``.

Dispatches to the compiled ``schubmult_cpp`` kernel when available, falling
back to ``_schubmult_q_fast_python`` (the ``medium_theta``-based recursion
with merged equal-length layers) otherwise.

**Arguments**:

- `perm_dict` - Mapping ``{Permutation: coeff}``.
- `v` - Permutation (or array-form list) indexing the quantum Schubert
  polynomial to multiply by.
- `q_var` - Generating set for the quantum parameters.
  

**Returns**:

- `dict` - Coefficient dict ``{Permutation: coeff}``, polynomial in ``q_var``.

<a id="schubmult.mult.quantum.schubmult_q"></a>

#### schubmult\_q

```python
def schubmult_q(perm_dict, v)
```

Multiply ``sum_u coeff_u S_u(x)`` by the quantum Schubert polynomial ``S_v``.

Reference (non-"fast") implementation: uses ``strict_theta`` and processes
every layer individually (no merging of equal-length adjacent layers), so it
is simpler but slower than ``schubmult_q_fast``. Results agree with
``schubmult_q_fast`` for all inputs.

**Arguments**:

- `perm_dict` - Mapping ``{Permutation: coeff}``.
- `v` - Permutation (or array-form list) indexing the quantum Schubert
  polynomial to multiply by.
  

**Returns**:

- `dict` - Coefficient dict ``{Permutation: coeff}``, polynomial in the
  default quantum parameters ``q``.

