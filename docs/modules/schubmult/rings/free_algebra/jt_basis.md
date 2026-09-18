<a id="schubmult.rings.free_algebra.jt_basis"></a>

# schubmult.rings.free\_algebra.jt\_basis

`JTBasis`: `JBasis` with a formal parameter ``t`` recording the number of stripped zeros.
Keys are ``(composition, power_of_t)``.

<a id="schubmult.rings.free_algebra.jt_basis.JTBasis"></a>

## JTBasis Objects

```python
class JTBasis(FreeAlgebraBasis)
```

JT basis of the free algebra (J basis with a parameter *t*).

Keys are ``(tuple, int)`` pairs where the tuple is a nonzero code
and the integer tracks a power of the parameter *t*.

<a id="schubmult.rings.free_algebra.jt_basis.JTBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.jt_basis.JTBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.jt_basis.JTBasis.from_perm"></a>

#### from\_perm

```python
@staticmethod
def from_perm(perm, n)
```

Extract a JT key from *perm* if the first *n* code entries are nonzero.

<a id="schubmult.rings.free_algebra.jt_basis.JTBasis.pare_schubert"></a>

#### pare\_schubert

```python
@staticmethod
def pare_schubert(perm)
```

Extract a nonzero trimcode from *perm*, or None if it contains zeros.

<a id="schubmult.rings.free_algebra.jt_basis.JTBasis.normalize_dct"></a>

#### normalize\_dct

```python
@staticmethod
def normalize_dct(dct)
```

Normalize a word dict by collecting zeros into leading positions.

<a id="schubmult.rings.free_algebra.jt_basis.JTBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return a *t*-weighted ``JT``-labelled display object.

<a id="schubmult.rings.free_algebra.jt_basis.JTBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from JTBasis to *other_basis*.

