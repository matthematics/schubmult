<a id="schubmult.symbolic.symmetric_polynomials.complete_sym"></a>

# schubmult.symbolic.symmetric\_polynomials.complete\_sym

Unevaluated (factorial) complete homogeneous symmetric polynomials as SymPy function atoms.

``H(p, k, xvars, yvars)`` (alias `FactorialCompleteSym`) is the factorial complete symmetric
polynomial of degree ``p`` in ``k`` generators with ``p + k - 1`` coefficient variables;
``h(p, k, xvars)`` (alias `CompleteSym`) is the ordinary ``h_p(x_1..x_k)``. The two families
are related by the duality ``H(p, k; x, y) = (-1)^p E(p, k + 1 - p; y, x)`` (`H.to_elem_sym`,
`H.from_elem_sym`), and divided differences are computed by passing through `E`.

<a id="schubmult.symbolic.symmetric_polynomials.complete_sym.CompleteSym_base"></a>

## CompleteSym\_base Objects

```python
class CompleteSym_base(Function)
```

Common behavior for `H` and `h`; ``expand_func`` evaluates via `complete_sym_poly`.

<a id="schubmult.symbolic.symmetric_polynomials.complete_sym.H"></a>

## H Objects

```python
class H(CompleteSym_base)
```

Factorial complete symmetric polynomial ``H(p, k, xvars, yvars)``; see the module docstring.

<a id="schubmult.symbolic.symmetric_polynomials.complete_sym.H.to_elem_sym"></a>

#### to\_elem\_sym

```python
def to_elem_sym()
```

``(-1)^p E(p, k + 1 - p; yvars, xvars)``: the same polynomial as a factorial elementary symmetric atom.

<a id="schubmult.symbolic.symmetric_polynomials.complete_sym.H.from_elem_sym"></a>

#### from\_elem\_sym

```python
@classmethod
def from_elem_sym(cls, elem, sign=False)
```

Inverse of `to_elem_sym`: the ``H`` atom equal to the `E` atom ``elem`` (with the ``(-1)^p`` if ``sign``).

<a id="schubmult.symbolic.symmetric_polynomials.complete_sym.H.split_out_vars"></a>

#### split\_out\_vars

```python
def split_out_vars(vars1, vars2=None)
```

``H_p(all) = sum_i H_i(vars1) H_{p-i}(rest)`` with the coefficient variables split accordingly.

<a id="schubmult.symbolic.symmetric_polynomials.complete_sym.H.divide_out_diff"></a>

#### divide\_out\_diff

```python
def divide_out_diff(v1, v2)
```

`E.divide_out_diff` transported through `to_elem_sym`/`from_elem_sym`.

<a id="schubmult.symbolic.symmetric_polynomials.complete_sym.H.from_expr_elem_sym"></a>

#### from\_expr\_elem\_sym

```python
@staticmethod
def from_expr_elem_sym(expr)
```

Replace every `E` atom in ``expr`` by the equal `H` atom.

<a id="schubmult.symbolic.symmetric_polynomials.complete_sym.H.to_expr_elem_sym"></a>

#### to\_expr\_elem\_sym

```python
@staticmethod
def to_expr_elem_sym(expr)
```

Replace every `H` atom in ``expr`` by the equal `E` atom.

<a id="schubmult.symbolic.symmetric_polynomials.complete_sym.H.div_diff"></a>

#### div\_diff

```python
def div_diff(v1, v2)
```

Divided difference, computed on the `E` side and converted back.

<a id="schubmult.symbolic.symmetric_polynomials.complete_sym.h"></a>

## h Objects

```python
class h(CompleteSym_base)
```

Ordinary complete homogeneous symmetric polynomial ``h(p, k, xvars)``.

