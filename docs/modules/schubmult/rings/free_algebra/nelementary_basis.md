<a id="schubmult.rings.free_algebra.nelementary_basis"></a>

# schubmult.rings.free\_algebra.nelementary\_basis

`NElementaryBasis`: the noncommutative elementary basis ``L`` of NSym inside the free algebra.

Keys are compositions (positive integers). ``L_alpha`` expands in words as the signed sum
``sum_{beta refines alpha} (-1)^(|alpha| - len(beta)) beta`` (refinements via SageMath), and the
product is concatenation.

<a id="schubmult.rings.free_algebra.nelementary_basis.NElementaryBasis"></a>

## NElementaryBasis Objects

```python
class NElementaryBasis(FreeAlgebraBasis)
```

Non-commutative elementary basis (L basis) of the free algebra.

Keys are tuples of positive integers. The transition to the word
basis uses composition refinements from SageMath.

<a id="schubmult.rings.free_algebra.nelementary_basis.NElementaryBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.nelementary_basis.NElementaryBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.nelementary_basis.NElementaryBasis.product"></a>

#### product

```python
@classmethod
def product(cls, key1, key2, coeff=S.One)
```

Concatenate two keys.

<a id="schubmult.rings.free_algebra.nelementary_basis.NElementaryBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return an ``L``-labelled display object for key *k*.

<a id="schubmult.rings.free_algebra.nelementary_basis.NElementaryBasis.transition_word"></a>

#### transition\_word

```python
@classmethod
def transition_word(cls, tup)
```

Transition an NElementary key to the word basis via composition refinements (requires SageMath).

<a id="schubmult.rings.free_algebra.nelementary_basis.NElementaryBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from NElementaryBasis to *other_basis*.

