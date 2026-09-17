<a id="schubmult.rings.free_algebra.fundamental_slide_basis"></a>

# schubmult.rings.free\_algebra.fundamental\_slide\_basis

<a id="schubmult.rings.free_algebra.fundamental_slide_basis.FundamentalSlideBasis"></a>

## FundamentalSlideBasis Objects

```python
class FundamentalSlideBasis(FreeAlgebraBasis)
```

Fundamental slide basis of the free algebra.

Keys are weak composition tuples. Transitions are computed via
polynomial algebra slide polynomial expansions.

<a id="schubmult.rings.free_algebra.fundamental_slide_basis.FundamentalSlideBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.fundamental_slide_basis.FundamentalSlideBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.fundamental_slide_basis.FundamentalSlideBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return an ``FS``-labelled display object for key *k*.

<a id="schubmult.rings.free_algebra.fundamental_slide_basis.FundamentalSlideBasis.transition_schubert"></a>

#### transition\_schubert

```python
@classmethod
def transition_schubert(cls, key)
```

Transition a fundamental slide key to the Schubert basis.

<a id="schubmult.rings.free_algebra.fundamental_slide_basis.FundamentalSlideBasis.transition_word"></a>

#### transition\_word

```python
@classmethod
@cache
def transition_word(cls, key)
```

Transition a fundamental slide key to the word basis.

<a id="schubmult.rings.free_algebra.fundamental_slide_basis.FundamentalSlideBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the FundamentalSlidePolyBasis as the dual.

<a id="schubmult.rings.free_algebra.fundamental_slide_basis.FundamentalSlideBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from FundamentalSlideBasis to *other_basis*.

