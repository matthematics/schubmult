<a id="schubmult.rings.free_algebra.z_basis"></a>

# schubmult.rings.free\_algebra.z\_basis

`ZBasis`: free-algebra basis indexed by compositions with no zeros, related to `SchubertBasis`
by shifting code entries by one and dropping zeros.

<a id="schubmult.rings.free_algebra.z_basis.ZBasis"></a>

## ZBasis Objects

```python
class ZBasis(FreeAlgebraBasis)
```

Z basis of the free algebra.

Keys are tuples of positive integers (no zeros). The Z basis is
related to the Schubert basis by incrementing/decrementing code
entries by 1 and dropping zeros.

<a id="schubmult.rings.free_algebra.z_basis.ZBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.z_basis.ZBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.z_basis.ZBasis.from_perm"></a>

#### from\_perm

```python
@staticmethod
def from_perm(perm, n)
```

Extract a Z basis key from *perm* if the first *n* code entries are nonzero.

<a id="schubmult.rings.free_algebra.z_basis.ZBasis.pare_schubert"></a>

#### pare\_schubert

```python
@staticmethod
def pare_schubert(perm)
```

Extract the nonzero trimcode of *perm*, or None if it contains interior zeros.

<a id="schubmult.rings.free_algebra.z_basis.ZBasis.product"></a>

#### product

```python
@classmethod
def product(cls, key1, key2, coeff=S.One)
```

Multiply two Z basis keys via shifted Schubert multiplication.

<a id="schubmult.rings.free_algebra.z_basis.ZBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return a ``Z``-labelled display object for key *k*.

<a id="schubmult.rings.free_algebra.z_basis.ZBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from ZBasis to *other_basis*.

