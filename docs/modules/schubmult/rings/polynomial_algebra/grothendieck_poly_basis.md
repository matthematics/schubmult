<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.grothendieck\_poly\_basis

`GrothendieckPolyBasis`: the Grothendieck polynomial basis of `PolynomialAlgebra`, at the
specialization ``beta = 1`` (without loss of generality: ``beta`` is recovered from the grading,
since the degree ``inv(w) + d`` part of ``G_w`` carries ``beta^d``).

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis"></a>

## GrothendieckPolyBasis Objects

```python
class GrothendieckPolyBasis(PolynomialBasis)
```

Grothendieck polynomial basis at ``beta = 1``.

Keys are ``(Permutation, length)`` pairs. Grothendieck polynomials form
the canonical basis for the polynomial algebra in Grothendieck calculus,
dual to the :class:`GrothendieckBasis` of the free algebra. The ``beta`` parameter
is set to 1 without loss of generality (see the module docstring).

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.product"></a>

#### product

```python
def product(key1, key2, coeff=S.One)
```

Multiply two Grothendieck keys using the Grothendieck ring multiplication.

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.transition_schubert"></a>

#### transition\_schubert

```python
def transition_schubert(dct)
```

Transition from Grothendieck basis to separated descents basis.

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.transition_glide_key"></a>

#### transition\_glide\_key

```python
def transition_glide_key(key)
```

Decompose a Grothendieck polynomial into glide polynomials via omega insertion on RC-graphs.

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.transition_glide"></a>

#### transition\_glide

```python
def transition_glide(dct)
```

Transition a Grothendieck dict to the glide polynomial basis.

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.to_monoms"></a>

#### to\_monoms

```python
def to_monoms(key)
```

Expand a Grothendieck key into a dict of monomial exponent tuples.

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual free algebra basis class (:class:`GrothendieckBasis`).

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.transition_grove_key"></a>

#### transition\_grove\_key

```python
def transition_grove_key(key)
```

Decompose a Grothendieck polynomial into grove polynomials via omega insertion on RC-graphs.

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.transition_grove"></a>

#### transition\_grove

```python
def transition_grove(dct)
```

Transition a Grothendieck dict to the grove polynomial basis.

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.transition_lascoux_key"></a>

#### transition\_lascoux\_key

```python
def transition_lascoux_key(key)
```

Decompose a Grothendieck polynomial into Lascoux polynomials via omega insertion on RC-graphs.

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.transition_lascoux"></a>

#### transition\_lascoux

```python
def transition_lascoux(dct)
```

Transition a Grothendieck dict to the Lascoux polynomial basis.

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from Grothendieck basis to *other_basis*.

