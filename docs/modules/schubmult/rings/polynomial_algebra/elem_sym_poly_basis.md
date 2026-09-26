<a id="schubmult.rings.polynomial_algebra.elem_sym_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.elem\_sym\_poly\_basis

`ElemSymPolyBasis`: the basis of products of elementary symmetric polynomials ``e_p(x_1..x_k)`` for `PolynomialAlgebra`.

<a id="schubmult.rings.polynomial_algebra.elem_sym_poly_basis.ElemSymPolyBasis"></a>

## ElemSymPolyBasis Objects

```python
class ElemSymPolyBasis(PolynomialBasis)
```

Elementary symmetric polynomial basis.

Keys are tuples encoding products of elementary symmetric polynomials
e_k(x_1, ..., x_n). Each key specifies degrees and variable counts
for the elementary symmetric factors.

<a id="schubmult.rings.polynomial_algebra.elem_sym_poly_basis.ElemSymPolyBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a ``(tuple/list, int)`` pair.

<a id="schubmult.rings.polynomial_algebra.elem_sym_poly_basis.ElemSymPolyBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* into a ``(tuple, int)`` key.

<a id="schubmult.rings.polynomial_algebra.elem_sym_poly_basis.ElemSymPolyBasis.transition_schubert"></a>

#### transition\_schubert

```python
def transition_schubert(dct)
```

Transition from elementary symmetric basis to Schubert basis.

<a id="schubmult.rings.polynomial_algebra.elem_sym_poly_basis.ElemSymPolyBasis.transition_monomial"></a>

#### transition\_monomial

```python
def transition_monomial(dct)
```

Transition from elementary symmetric basis to monomial basis.

<a id="schubmult.rings.polynomial_algebra.elem_sym_poly_basis.ElemSymPolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from this basis to *other_basis*.

