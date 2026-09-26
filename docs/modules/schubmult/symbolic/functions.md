<a id="schubmult.symbolic.functions"></a>

# schubmult.symbolic.functions

SymEngine-first wrappers (`expand`, `symbols`, `sympify`) that fall back to SymPy, plus small helpers.

<a id="schubmult.symbolic.functions.expand"></a>

#### expand

```python
def expand(obj, **kwargs)
```

Expand with SymEngine; use SymPy if keyword options are given or SymEngine fails.

<a id="schubmult.symbolic.functions.symbols"></a>

#### symbols

```python
def symbols(*args, **kwargs)
```

SymEngine ``symbols``.

<a id="schubmult.symbolic.functions.sympify"></a>

#### sympify

```python
def sympify(val)
```

SymEngine ``sympify``, falling back to SymPy for objects SymEngine cannot convert.

<a id="schubmult.symbolic.functions.is_of_func_type"></a>

#### is\_of\_func\_type

```python
def is_of_func_type(elem, typ)
```

``isinstance`` that also sees through SymEngine ``PyFunction`` wrappers around SymPy functions.

<a id="schubmult.symbolic.functions.expand_seq"></a>

#### expand\_seq

```python
def expand_seq(seq, genset)
```

The monomial ``genset[1]**seq[0] * genset[2]**seq[1] * ...`` (1-indexed generators).

<a id="schubmult.symbolic.functions.prod"></a>

#### prod

```python
def prod(a, start=1)
```

Product of the elements of ``a`` times ``start`` (same as ``sympy.prod``).

<a id="schubmult.symbolic.functions.efficient_subs"></a>

#### efficient\_subs

```python
def efficient_subs(expr, subs_dict)
```

``expr.subs`` restricted to the entries of ``subs_dict`` that actually occur in ``expr``.

