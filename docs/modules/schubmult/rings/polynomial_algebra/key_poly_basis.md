<a id="schubmult.rings.polynomial_algebra.key_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.key\_poly\_basis

`KeyPolyBasis`: the key polynomial (Demazure character) basis of `PolynomialAlgebra`, indexed by weak compositions.

<a id="schubmult.rings.polynomial_algebra.key_poly_basis.traverse_demaz"></a>

#### traverse\_demaz

```python
def traverse_demaz(pl, w)
```

Traverse the Demazure graph starting from plactic element *pl* along the code word of *w*.

Yields all distinct elements reachable by successive lowering operators.

<a id="schubmult.rings.polynomial_algebra.key_poly_basis.KeyPolyBasis"></a>

## KeyPolyBasis Objects

```python
class KeyPolyBasis(PolynomialBasis)
```

Key polynomial (Demazure character) basis.

Keys are weak compositions. Key polynomials are characters of Demazure
modules, computed by applying Demazure operators to highest-weight
plactic tableaux.

<a id="schubmult.rings.polynomial_algebra.key_poly_basis.KeyPolyBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual free algebra basis class (:class:`KeyBasis`).

<a id="schubmult.rings.polynomial_algebra.key_poly_basis.KeyPolyBasis.to_monoms"></a>

#### to\_monoms

```python
def to_monoms(key)
```

Expand a key polynomial key into a dict of monomial exponent tuples.

<a id="schubmult.rings.polynomial_algebra.key_poly_basis.KeyPolyBasis.expand"></a>

#### expand

```python
def expand(dct)
```

Expand a key basis dict into a symbolic polynomial expression.

<a id="schubmult.rings.polynomial_algebra.key_poly_basis.KeyPolyBasis.transition_monomial"></a>

#### transition\_monomial

```python
def transition_monomial(dct)
```

Transition from key basis to monomial basis.

<a id="schubmult.rings.polynomial_algebra.key_poly_basis.KeyPolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from key basis to *other_basis*.

