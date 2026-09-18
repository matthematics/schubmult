<a id="schubmult.symbolic.poly.variables"></a>

# schubmult.symbolic.poly.variables

Generating sets: indexed families of variables ``x_0, x_1, x_2, ...`` used as ring generators.

`GeneratingSet("x")` interns ``DEF_GENSET_SIZE`` symbols ``x_0..x_99``; ``gs[i]`` is the symbol
``x_i``, and the polynomial variables are ``x_1, x_2, ...`` (``x_0`` is unused), so an exponent
tuple ``(a_1, ..., a_n)`` means ``x_1^{a_1} ... x_n^{a_n}``. `MaskedGeneratingSet` hides a set of indices of a base set,
`CustomGeneratingSet` wraps an arbitrary sequence of expressions, and `ZeroGeneratingSet`
returns 0 for every index (used for single Schubert polynomials as a degenerate coefficient
set). `genset_dict_from_expr` converts a polynomial expression into ``{exponent_tuple: coeff}``.

<a id="schubmult.symbolic.poly.variables.GeneratingSet_base"></a>

## GeneratingSet\_base Objects

```python
class GeneratingSet_base()
```

Interface for generating sets: indexing, length, ``index(symbol)`` (``-1`` if absent), and ``label``.

<a id="schubmult.symbolic.poly.variables.ZeroGeneratingSet"></a>

## ZeroGeneratingSet Objects

```python
class ZeroGeneratingSet(GeneratingSet_base)
```

A generating set every entry of which is ``0``; contains no symbols.

<a id="schubmult.symbolic.poly.variables.GeneratingSet"></a>

## GeneratingSet Objects

```python
class GeneratingSet(GeneratingSet_base)
```

The interned family ``name_0, name_1, ...``; ``gs[i]`` is the symbol ``name_i`` and ``gs(i)`` is ``gs[i - 1]``.

<a id="schubmult.symbolic.poly.variables.GeneratingSet.__call__"></a>

#### \_\_call\_\_

```python
def __call__(index)
```

1-indexed

<a id="schubmult.symbolic.poly.variables.GeneratingSet.label"></a>

#### label

```python
@property
def label()
```

The variable name, e.g. ``"x"``.

<a id="schubmult.symbolic.poly.variables.GeneratingSet.index"></a>

#### index

```python
def index(v)
```

Position of the symbol ``v`` in this set, or ``-1``.

<a id="schubmult.symbolic.poly.variables.MaskedGeneratingSet"></a>

## MaskedGeneratingSet Objects

```python
class MaskedGeneratingSet(GeneratingSet_base)
```

A base generating set with the (1-indexed) positions in ``index_mask`` removed and the rest
renumbered consecutively; ``complement()`` gives the set of the masked variables instead.

<a id="schubmult.symbolic.poly.variables.MaskedGeneratingSet.base_genset"></a>

#### base\_genset

```python
@property
def base_genset()
```

The underlying unmasked generating set.

<a id="schubmult.symbolic.poly.variables.MaskedGeneratingSet.index_mask"></a>

#### index\_mask

```python
@property
def index_mask()
```

Sorted tuple of the hidden 1-indexed positions.

<a id="schubmult.symbolic.poly.variables.MaskedGeneratingSet.complement"></a>

#### complement

```python
def complement()
```

The masked set on the complementary positions.

<a id="schubmult.symbolic.poly.variables.MaskedGeneratingSet.__call__"></a>

#### \_\_call\_\_

```python
def __call__(index)
```

1-indexed

<a id="schubmult.symbolic.poly.variables.CustomGeneratingSet"></a>

## CustomGeneratingSet Objects

```python
class CustomGeneratingSet(GeneratingSet_base)
```

A generating set over an explicit sequence of expressions (sympified on construction).

<a id="schubmult.symbolic.poly.variables.CustomGeneratingSet.__call__"></a>

#### \_\_call\_\_

```python
def __call__(index)
```

1-indexed

<a id="schubmult.symbolic.poly.variables.NotEnoughGeneratorsError"></a>

## NotEnoughGeneratorsError Objects

```python
class NotEnoughGeneratorsError(ValueError)
```

Raised when an operation needs more generators than a generating set provides.

<a id="schubmult.symbolic.poly.variables.poly_genset"></a>

#### poly\_genset

```python
@cache
def poly_genset(v: str)
```

``GeneratingSet(v)``, or a `ZeroGeneratingSet` for the sentinels ``ZeroVar``/``NoneVar``.

<a id="schubmult.symbolic.poly.variables.genset_dict_from_expr"></a>

#### genset\_dict\_from\_expr

```python
def genset_dict_from_expr(expr, genset, length=None)
```

Write a polynomial in the generators of ``genset`` as ``{exponent_tuple: coeff}``.

Exponent tuples are 0-indexed by generator position ``genset(i) -> tuple[i - 1]`` and have
length ``length`` (default: the largest generator index present). Factors free of the
generators go into the coefficient; a factor mixing generators with other symbols raises.

