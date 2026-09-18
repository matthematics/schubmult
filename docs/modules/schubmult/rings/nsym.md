<a id="schubmult.rings.nsym"></a>

# schubmult.rings.nsym

`NSym`: a `FreeAlgebra` indexed by compositions, multiplied through the separated-descents Schubert ring.

A key is a composition ``alpha`` (a tuple of positive integers), printed ``N(alpha)``, and is
identified with the Schubert key ``(uncode(alpha - 1), len(alpha))`` of
`schubmult.rings.schubert.separated_descents.SeparatedDescentsRing` via `NSym.sepify` /
`NSym.from_sep`. The product is the separated-descents product transported back to
compositions; when ``FreeAlgebra.CAP`` is set the result is truncated to keys of at most that
length. Right multiplication by a Schubert element acts by the skew operation ``/``.

<a id="schubmult.rings.nsym.NSym"></a>

## NSym Objects

```python
class NSym(FreeAlgebra)
```

Free algebra on compositions with the separated-descents product. See the module docstring.

<a id="schubmult.rings.nsym.NSym.__init__"></a>

#### \_\_init\_\_

```python
def __init__(domain=None)
```

Create the ring over ``domain`` (default ``EXRAW``); the empty composition is the identity.

<a id="schubmult.rings.nsym.NSym.printing_term"></a>

#### printing\_term

```python
def printing_term(k)
```

Display as ``N(alpha)``.

<a id="schubmult.rings.nsym.NSym.rmul"></a>

#### rmul

```python
def rmul(elem, other)
```

Scale coefficients by the scalar ``other``.

<a id="schubmult.rings.nsym.NSym.sepify"></a>

#### sepify

```python
def sepify(elem)
```

Map ``alpha -> (uncode(alpha - 1), len(alpha))`` into the separated-descents Schubert ring.

<a id="schubmult.rings.nsym.NSym.from_sep"></a>

#### from\_sep

```python
def from_sep(elem)
```

Inverse of `sepify`: pad or cut the code of ``perm`` to length ``n`` and add 1 to each entry.

<a id="schubmult.rings.nsym.NSym.mul"></a>

#### mul

```python
def mul(elem, other)
```

Scalar multiplication, or the separated-descents product of two elements (truncated by
``FreeAlgebra.CAP`` if set).

<a id="schubmult.rings.nsym.NSymElement"></a>

## NSymElement Objects

```python
class NSymElement(FreeAlgebraElement)
```

Element of `NSym`: a dict from compositions to coefficients with SymPy-compatible printing.

<a id="schubmult.rings.nsym.NSymElement.__rmul__"></a>

#### \_\_rmul\_\_

```python
def __rmul__(other)
```

Scalar on the left, or a Schubert element acting by the skew operation ``self / perm``.

