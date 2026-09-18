<a id="schubmult.symbolic.symmetric_polynomials.elem_sym"></a>

# schubmult.symbolic.symmetric\_polynomials.elem\_sym

Unevaluated (factorial) elementary symmetric polynomials as SymPy function atoms.

``E(p, k, xvars, yvars)`` (alias `FactorialElemSym`) is the factorial elementary symmetric
polynomial of degree ``p`` in the ``k`` generators ``xvars`` with coefficient variables
``yvars`` (``k + 1 - p`` of them are used); ``e(p, k, xvars)`` (alias `ElemSym`) is the ordinary
``e_p(x_1..x_k)``. Both stay symbolic so Schubert polynomials can be manipulated in the SEM/CEM
bases; ``expand_func`` evaluates them via `schubmult.symbolic.poly.schub_poly.elem_sym_poly`.
They implement divided differences (`div_diff`, `divide_out_diff`), variable splitting
(`split_out_vars`), and canonicalize on construction (``E(p, k, ...) = 0`` if ``p > k``, ``1`` if
``p == 0``, and a shared variable between the two sets cancels).

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym.ElemSym_base"></a>

## ElemSym\_base Objects

```python
class ElemSym_base(Function)
```

Common behavior for `E` and `e`: substitution acts only on the variable arguments, and the
``degree``/``numvars``/``genvars``/``coeffvars`` accessors expose the parameters.

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym.ElemSym_base.degree"></a>

#### degree

```python
@property
def degree()
```

``p``.

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym.ElemSym_base.numvars"></a>

#### numvars

```python
@property
def numvars()
```

``k``, the number of generators.

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym.ElemSym_base.genvars"></a>

#### genvars

```python
@property
def genvars()
```

The ``x`` variables (sorted tuple).

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym.ElemSym_base.coeffvars"></a>

#### coeffvars

```python
@property
def coeffvars()
```

The ``y`` (coefficient) variables.

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym.E"></a>

## E Objects

```python
class E(ElemSym_base)
```

Factorial elementary symmetric polynomial ``E(p, k, xvars, yvars)``; see the module docstring.

Variables may be passed as two iterables or flattened (``k`` x's followed by ``k + 1 - p`` y's).

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym.E.cauchy"></a>

#### cauchy

```python
@staticmethod
def cauchy(fnc, genset)
```

Rewrite ``fnc`` so its coefficient variables are the initial segment of ``genset``, using
``E(p,k;..y..) = E(p,k;..y'..) + (y - y') E(p-1,k-1;..y'..)`` one variable at a time.

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym.E.split_out_vars"></a>

#### split\_out\_vars

```python
def split_out_vars(vars1, vars2=None)
```

Split the generators into ``vars1`` and the rest: ``E(p, k) = sum_i E(i, k1; ..) E(p - i, k2; ..)``
with the coefficient variables distributed accordingly. With ``vars1=None`` the split is
made on the coefficient variables ``vars2`` instead.

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym.E.divide_out_diff"></a>

#### divide\_out\_diff

```python
def divide_out_diff(v1, v2)
```

``(self - self|_{v1 -> v2}) / (v1 - v2)`` in closed form: removing a generator lowers ``p`` and
``k`` by one; a coefficient variable gives the corresponding signed term.

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym.E.div_diff"></a>

#### div\_diff

```python
def div_diff(v1, v2)
```

Divided difference ``partial_{v1, v2}`` in closed form (antisymmetric in ``v1``, ``v2``).

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym.E.pull_out_vars"></a>

#### pull\_out\_vars

```python
def pull_out_vars(var1, var2, min_degree=1)
```

Write ``self = self|_{var1 -> var2} + (var1 - var2) * divide_out_diff(var1, var2)`` when ``var1``
is a generator and ``var2`` a coefficient variable (and ``p >= min_degree``).

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym.e"></a>

## e Objects

```python
class e(ElemSym_base)
```

Ordinary elementary symmetric polynomial ``e(p, k, xvars)``; see the module docstring.

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym.e.split_out_vars"></a>

#### split\_out\_vars

```python
def split_out_vars(vars1, vars2=None)
```

``e_p(all) = sum_i e_i(vars1) e_{p-i}(rest)``.

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym.e.coeffvars"></a>

#### coeffvars

```python
@property
def coeffvars()
```

No coefficient variables: a `ZeroGeneratingSet`.

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym.e.divide_out_diff"></a>

#### divide\_out\_diff

```python
def divide_out_diff(v1, v2)
```

``e_{p-1}`` of the remaining generators if ``v1`` is a generator, else 0.

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym.e.div_diff"></a>

#### div\_diff

```python
def div_diff(v1, v2)
```

Divided difference: ``+/- e_{p-1}`` of the remaining generators, or 0 if neither variable is a generator.

