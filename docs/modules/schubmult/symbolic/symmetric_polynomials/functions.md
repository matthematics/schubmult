<a id="schubmult.symbolic.symmetric_polynomials.functions"></a>

# schubmult.symbolic.symmetric\_polynomials.functions

Expression-level utilities for symbolic elementary symmetric polynomials.

The accessors `genvars`/`coeffvars`/`degree`/`numvars` see through SymEngine ``PyFunction``
wrappers; `split_out_vars`/`pull_out_vars` map the corresponding `E` methods over a whole
expression tree; `canonicalize_elem_syms` rewrites products of ``E`` factors into a normal form
(each factor of full degree ``p == k``, grouped by first coefficient variable).

<a id="schubmult.symbolic.symmetric_polynomials.functions.genvars"></a>

#### genvars

```python
def genvars(obj)
```

``obj.genvars``, unwrapping a SymEngine ``PyFunction`` if needed.

<a id="schubmult.symbolic.symmetric_polynomials.functions.coeffvars"></a>

#### coeffvars

```python
def coeffvars(obj)
```

``obj.coeffvars``, unwrapping a SymEngine ``PyFunction`` if needed.

<a id="schubmult.symbolic.symmetric_polynomials.functions.degree"></a>

#### degree

```python
def degree(obj)
```

The degree ``p`` of an elementary symmetric atom (unwrapping if needed).

<a id="schubmult.symbolic.symmetric_polynomials.functions.numvars"></a>

#### numvars

```python
def numvars(obj)
```

The variable count ``k`` of an elementary symmetric atom (unwrapping if needed).

<a id="schubmult.symbolic.symmetric_polynomials.functions.canonicalize_elem_syms"></a>

#### canonicalize\_elem\_syms

```python
def canonicalize_elem_syms(expr, combine_equal=False)
```

Normal form for expressions in `FactorialElemSym`: split every factor with ``p < k`` in half
until all factors have ``p == k``, then within each product regroup factors sharing a first
coefficient variable (merging them into one factor if ``combine_equal``).

<a id="schubmult.symbolic.symmetric_polynomials.functions.canonicalize_elem_syms_coeff"></a>

#### canonicalize\_elem\_syms\_coeff

```python
def canonicalize_elem_syms_coeff(expr, combine_equal=False)
```

`canonicalize_elem_syms` splitting on coefficient variables instead of generators.

<a id="schubmult.symbolic.symmetric_polynomials.functions.split_out_vars"></a>

#### split\_out\_vars

```python
def split_out_vars(expr, vars1, vars2)
```

Apply ``split_out_vars(vars1, vars2)`` to every elementary symmetric atom in ``expr``.

<a id="schubmult.symbolic.symmetric_polynomials.functions.pull_out_vars"></a>

#### pull\_out\_vars

```python
def pull_out_vars(expr, var1, var2, min_degree=1)
```

Apply ``pull_out_vars(var1, var2, min_degree)`` to every elementary symmetric atom in ``expr``.

<a id="schubmult.symbolic.symmetric_polynomials.functions.elem_sym_unify"></a>

#### elem\_sym\_unify

```python
def elem_sym_unify(expr, arg=None)
```

Recursively walk ``expr`` (currently a structural no-op; the pattern-based unification is
commented out).

