<a id="schubmult.rings.free_algebra.monomial_slide_basis"></a>

# schubmult.rings.free\_algebra.monomial\_slide\_basis

`MonomialSlideBasis`: the free-algebra basis dual to monomial slide polynomials. Keys are
weak compositions; transitions use coarsenings of compositions.

<a id="schubmult.rings.free_algebra.monomial_slide_basis.MonomialSlideBasis"></a>

## MonomialSlideBasis Objects

```python
class MonomialSlideBasis(FreeAlgebraBasis)
```

Monomial slide basis of the free algebra.

Keys are weak composition tuples. Transitions use monomial slide
polynomial expansions and coarsenings of compositions.

<a id="schubmult.rings.free_algebra.monomial_slide_basis.MonomialSlideBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.monomial_slide_basis.MonomialSlideBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.monomial_slide_basis.MonomialSlideBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return an ``MS``-labelled display object for key *k*.

<a id="schubmult.rings.free_algebra.monomial_slide_basis.MonomialSlideBasis.transition_fundamental_slide"></a>

#### transition\_fundamental\_slide

```python
@classmethod
@cache
def transition_fundamental_slide(cls, key)
```

Transition a monomial slide key to the fundamental slide basis.

<a id="schubmult.rings.free_algebra.monomial_slide_basis.MonomialSlideBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from MonomialSlideBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.monomial_slide_basis.MonomialSlideBasis.transition_word"></a>

#### transition\_word

```python
@classmethod
@cache
def transition_word(cls, key)
```

Transition a monomial slide key to the word basis.

