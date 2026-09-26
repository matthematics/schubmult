<a id="schubmult.mult.quantum_double"></a>

# schubmult.mult.quantum\_double

Quantum double Schubert polynomial multiplication.

Implements ``schubmult_q_double``/``schubmult_q_double_fast``: the product of a
linear combination of quantum double Schubert polynomials ``S_u(x, var2)`` with
a single ``S_v(x, var3)``, returned as a coefficient dict ``{w: coeff}``
polynomial in ``var2``, ``var3``, and the quantum parameters ``q_1, q_2, ...``.
Uses the same ``theta``/v-path recursion as ``schubmult.mult.double``, with the
elementary-symmetric step generalized to the quantum moves of
``elem_sym_perms_q`` and the coefficient function to ``elem_sym_func_q``.

Also provides: ``mult_poly_q_double`` (multiply by an arbitrary polynomial),
``apply_peterson_woodward`` (parabolic quantum via the Peterson-Woodward
comparison theorem), ``q_posify``/``q_partial_posify_generic`` (manifestly
positive display of quantum structure constants), ``schubpoly_quantum``
(the quantum Schubert polynomial itself), ``nil_hecke`` (quantum nilHecke
action), and ``factor_out_q`` (split a polynomial by its ``q``-monomials).

<a id="schubmult.mult.quantum_double.single_variable"></a>

#### single\_variable

```python
def single_variable(coeff_dict, varnum, var_y=None, q_var=_vars.q_var)
```

Multiply ``sum_u coeff_u S_u(x, var_y)`` by the single variable ``x_varnum`` (quantum equivariant Monk rule).

The diagonal term contributes ``var_y[u(varnum)]``; the off-diagonal terms come from
``elem_sym_positional_perms_q``, each carrying its ``q``-monomial and sign.

<a id="schubmult.mult.quantum_double.mult_poly_q_double"></a>

#### mult\_poly\_q\_double

```python
def mult_poly_q_double(coeff_dict,
                       poly,
                       var_x=None,
                       var_y=None,
                       q_var=_vars.q_var)
```

Multiply ``sum_u coeff_u S_u(x, var_y)`` by an arbitrary polynomial ``poly`` in ``var_x``.

Recurses over the ``Add``/``Mul``/``Pow`` structure of ``poly``, dispatching
single-variable leaves to ``single_variable``; mirrors ``mult_poly_double``.

<a id="schubmult.mult.quantum_double.mult_poly_q_double_alt"></a>

#### mult\_poly\_q\_double\_alt

```python
def mult_poly_q_double_alt(coeff_dict,
                           poly,
                           var_x=None,
                           var_y=None,
                           q_var=_vars.q_var)
```

Variant of ``mult_poly_q_double`` that folds each factor via ``schubmult_q_double_dict_fast``
instead of ``single_variable``.

<a id="schubmult.mult.quantum_double.nil_hecke"></a>

#### nil\_hecke

```python
def nil_hecke(perm_dict, v, n, var2=None, var3=None)
```

Quantum nilHecke action: like ``schubmult_q_double`` but using the descent-side
``elem_sym_perms_q_op`` moves (bounded by ``n``) with ``up`` and ``up2`` swapped in the
coefficient function.

<a id="schubmult.mult.quantum_double.schubmult_q_double_pair"></a>

#### schubmult\_q\_double\_pair

```python
@cache
def schubmult_q_double_pair(perm1, perm2, var2=None, var3=None, q_var=None)
```

``schubmult_q_double_fast`` specialized to a single ``perm1`` with coefficient 1, cached.

<a id="schubmult.mult.quantum_double.schubmult_q_double_pair_generic"></a>

#### schubmult\_q\_double\_pair\_generic

```python
@cache
def schubmult_q_double_pair_generic(perm1, perm2)
```

``schubmult_q_double_pair`` with the fixed generic alphabets ``_vars.var_g1``/``_vars.var_g2``/``_vars.q_var``.

<a id="schubmult.mult.quantum_double.schubmult_q_generic_partial_posify"></a>

#### schubmult\_q\_generic\_partial\_posify

```python
@cache
def schubmult_q_generic_partial_posify(u2, v2)
```

Manifestly positive (where possible) expansion of ``S_{u2} * S_{v2}`` over the generic alphabets,
applying ``q_partial_posify_generic`` to each coefficient.

<a id="schubmult.mult.quantum_double.q_posify"></a>

#### q\_posify

```python
def q_posify(u, v, w, val, var2, var3, q_var, msg)
```

Manifestly positive representation of the quantum double structure constant ``c^w_{u,v}``.

Splits ``val`` by ``q``-monomial (``factor_out_q``), then for each piece either takes it
as-is (integer, or when ``v``'s inverse code is already in medium-theta form), reduces the
triple ``(u, v, w)`` via ``reduce_q_coeff`` until the ``q``-monomial becomes trivial and
delegates to the classical ``posify``, or falls back to ``compute_positive_rep``.
Raises if the reconstruction does not equal ``val``.

<a id="schubmult.mult.quantum_double.q_partial_posify_generic"></a>

#### q\_partial\_posify\_generic

```python
def q_partial_posify_generic(val, u, v, w)
```

Like ``q_posify`` over the generic alphabets, but only attempts positivity when ``v`` contains
a ``1432`` or ``312`` pattern (otherwise the raw value is already manifestly positive), and
leaves non-reducible ``q``-pieces unchanged rather than running the LP.

<a id="schubmult.mult.quantum_double.apply_peterson_woodward"></a>

#### apply\_peterson\_woodward

```python
def apply_peterson_woodward(coeff_dict,
                            parabolic_index,
                            q_var=_vars.q_var,
                            n=None)
```

Project a full-flag quantum product onto the parabolic quantum cohomology for ``parabolic_index``.

Implements the Peterson-Woodward comparison: for each ``q``-monomial of each coefficient,
checks the ``omega``/``check_blocks`` compatibility conditions on the exponent vector,
multiplies the indexing permutation by the appropriate parabolic longest elements, keeps
only the ``parabolic``-minimal results, and reindexes the surviving ``q`` variables.

**Arguments**:

- `coeff_dict` - Full-flag quantum coefficient dict ``{Permutation: coeff}``.
- `parabolic_index` - Sorted list of 1-indexed positions generating the parabolic subgroup.
- `q_var` - Quantum parameter generating set.
- `n` - Ambient flag size ``S_n``; results indexed by longer permutations are dropped. Defaults to
  ``parabolic_index[-1] + 1``, which undercounts when the last block has size 1 (it
  contributes no reflection), so callers that know the block sizes should pass their sum.
  

**Returns**:

- `dict` - Parabolic quantum coefficient dict ``{Permutation: coeff}``.

<a id="schubmult.mult.quantum_double.elem_sym_func_q_q"></a>

#### elem\_sym\_func\_q\_q

```python
def elem_sym_func_q_q(k,
                      i,
                      u1,
                      u2,
                      v1,
                      v2,
                      udiff,
                      vdiff,
                      varl1,
                      varl2,
                      q_var=_vars.q_var)
```

Fully-quantum coefficient function for the v-path recursion (used by ``schubpoly_quantum``):
the quantum elementary symmetric polynomial ``elem_sym_poly_q`` in the fixed-window ``y``
variables and the ``call_zvars`` ``z`` variables.

<a id="schubmult.mult.quantum_double.schubpoly_quantum"></a>

#### schubpoly\_quantum

```python
def schubpoly_quantum(v, var_x=None, var_y=None, q_var=_vars.q_var, coeff=1)
```

The quantum double Schubert polynomial ``S_v(var_x, var_y)`` itself, as a symbolic expression.

Runs the v-path recursion starting from the identity with ``elem_sym_func_q_q`` and reads off
the coefficient of the identity permutation.

<a id="schubmult.mult.quantum_double.schubmult_q_double"></a>

#### schubmult\_q\_double

```python
def schubmult_q_double(perm_dict, v, var2=None, var3=None, q_var=_vars.q_var)
```

Multiply ``sum_u coeff_u S_u(x, var2)`` by the quantum double Schubert polynomial ``S_v(x, var3)``.

Reference (non-"fast") implementation: uses ``strict_theta`` and processes every layer
individually. Results agree with ``schubmult_q_double_fast``.

**Arguments**:

- `perm_dict` - Mapping ``{Permutation: coeff}``.
- `v` - Permutation to multiply by.
- `var2` - Secondary alphabet attached to ``perm_dict``'s permutations.
- `var3` - Secondary alphabet attached to ``v``.
- `q_var` - Quantum parameter generating set.
  

**Returns**:

- `dict` - Coefficient dict ``{Permutation: coeff}``.

<a id="schubmult.mult.quantum_double.schubmult_q_double_dict_fast"></a>

#### schubmult\_q\_double\_dict\_fast

```python
def schubmult_q_double_dict_fast(perm_dict1,
                                 perm_dict2,
                                 var2=None,
                                 var3=None,
                                 q_var=_vars.q_var)
```

Multiply two coefficient dicts of quantum double Schubert polynomials together.

Sums ``schubmult_q_double_fast(perm_dict1, v, ...)`` scaled by ``coeff2_v`` over ``v`` in ``perm_dict2``.

<a id="schubmult.mult.quantum_double.schubmult_q_double_fast"></a>

#### schubmult\_q\_double\_fast

```python
def schubmult_q_double_fast(perm_dict,
                            v,
                            var2=None,
                            var3=None,
                            q_var=_vars.q_var)
```

Multiply ``sum_u coeff_u S_u(x, var2)`` by the quantum double Schubert polynomial ``S_v(x, var3)``.

Dispatches to the compiled ``schubmult_cpp`` kernel when available (and both secondary
alphabets are given), falling back to ``_schubmult_q_double_fast_python`` (the
``medium_theta``-based recursion with merged equal-length layers) otherwise.

**Arguments**:

- `perm_dict` - Mapping ``{Permutation: coeff}``.
- `v` - Permutation to multiply by.
- `var2` - Secondary alphabet attached to ``perm_dict``'s permutations.
- `var3` - Secondary alphabet attached to ``v``.
- `q_var` - Quantum parameter generating set.
  

**Returns**:

- `dict` - Coefficient dict ``{Permutation: coeff}``.

<a id="schubmult.mult.quantum_double.sum_q_dict"></a>

#### sum\_q\_dict

```python
def sum_q_dict(q_dict1, q_dict2)
```

Add two ``{q_monomial: coeff}`` dicts.

<a id="schubmult.mult.quantum_double.mul_q_dict"></a>

#### mul\_q\_dict

```python
def mul_q_dict(q_dict1, q_dict2)
```

Multiply two ``{q_monomial: coeff}`` dicts (convolution over monomials).

<a id="schubmult.mult.quantum_double.factor_out_q"></a>

#### factor\_out\_q

```python
def factor_out_q(poly, q_var=_vars.q_var)
```

Split ``poly`` by its ``q``-monomials: return ``{q_monomial: coefficient}`` with coefficients
free of ``q_var`` variables. Recurses over the ``Add``/``Mul``/``Pow`` structure; a polynomial
with no ``q`` variables maps to ``{1: poly}``.

