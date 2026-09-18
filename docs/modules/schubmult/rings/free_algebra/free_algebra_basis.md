<a id="schubmult.rings.free_algebra.free_algebra_basis"></a>

# schubmult.rings.free\_algebra.free\_algebra\_basis

`FreeAlgebraBasis`: the interface a basis must implement to plug into `FreeAlgebra`.

A basis is a *class* (its methods are classmethods) defining a key type (``is_key``/
``as_key``/``zero_monom``), how to print a key, ``transition(other_basis)`` returning a
key -> ``{key: coeff}`` function into another basis, and ``dual_basis()`` naming the
`schubmult.rings.polynomial_algebra` basis it is dual to. Products, coproducts, and
the word-level operations all have default implementations that route through the
`WordBasis` via ``compose_transition``.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis"></a>

## FreeAlgebraBasis Objects

```python
class FreeAlgebraBasis()
```

Abstract base for free-algebra bases; see the module docstring.

Subclasses override the key methods (``is_key``, ``as_key``, ``zero_monom``,
``printing_term``, ``transition``, ``dual_basis``) and may override ``product``/
``coproduct`` with a direct rule; otherwise everything is computed in the `WordBasis`
and transported back.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a valid key for this basis.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.from_rc_graph"></a>

#### from\_rc\_graph

```python
@classmethod
def from_rc_graph(cls, rc_graph)
```

Convert an RC graph to a basis-keyed dict.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* into a canonical key for this basis.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a function mapping keys of this basis to dicts in *other_basis*.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, key)
```

Return the display symbol for *key*.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.compose_transition"></a>

#### compose\_transition

```python
@classmethod
def compose_transition(cls, tkeyfunc, output)
```

Apply a key-level transition function to each key in *output*.

For each ``(key, v)`` in *output*, expands ``tkeyfunc(key)`` and
accumulates the results weighted by *v*.

**Arguments**:

- `tkeyfunc` - A function mapping a key to a ``{key: coeff}`` dict.
- `output` - A ``{key: coeff}`` dict to transform.
  

**Returns**:

  A merged ``{key: coeff}`` dict in the target basis.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

The `schubmult.rings.polynomial_algebra` basis this basis is dual to under the word/monomial pairing.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.change_tensor_basis"></a>

#### change\_tensor\_basis

```python
@classmethod
def change_tensor_basis(cls, tensor_elem, basis1, basis2)
```

Change the bases of both factors of a tensor element.

**Arguments**:

- `tensor_elem` - An element of a tensor product ring.
- `basis1` - Target basis for the left factor.
- `basis2` - Target basis for the right factor.
  

**Returns**:

  The tensor element re-expressed in the new bases.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.coproduct"></a>

#### coproduct

```python
@classmethod
@cache
def coproduct(cls, key)
```

Compute the coproduct of *key* by delegating through WordBasis.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.bcoproduct"></a>

#### bcoproduct

```python
@classmethod
@cache
def bcoproduct(cls, key)
```

Compute the bar-coproduct of *key* by delegating through WordBasis.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.product"></a>

#### product

```python
@classmethod
@cache
def product(cls, key1, key2, coeff=S.One)
```

Multiply two keys by transitioning to WordBasis and back.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.internal_product"></a>

#### internal\_product

```python
@classmethod
def internal_product(cls, key1, key2, coeff=S.One)
```

The internal (Kronecker) product of NSym (see `WordBasis.internal_product`), computed via the word basis.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.inject"></a>

#### inject

```python
@classmethod
def inject(cls, key1, i, key2, coeff=S.One)
```

Inject *key2* into *key1* at position *i* by delegating through WordBasis.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.prefix"></a>

#### prefix

```python
@classmethod
def prefix(cls, key, length, coeff=S.One)
```

Extract a prefix of *length* letters by delegating through WordBasis.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.suffix"></a>

#### suffix

```python
@classmethod
def suffix(cls, key, length, coeff=S.One)
```

Extract a suffix of *length* letters by delegating through WordBasis.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.interval"></a>

#### interval

```python
@classmethod
def interval(cls, key, start, stop, coeff=S.One)
```

Extract a subword from *start* to *stop* by delegating through WordBasis.

<a id="schubmult.rings.free_algebra.free_algebra_basis.__getattr__"></a>

#### \_\_getattr\_\_

```python
def __getattr__(name)
```

Lazily resolve basis classes to avoid circular imports between basis modules.

