<a id="schubmult.rings.free_algebra.grothendieck_basis"></a>

# schubmult.rings.free\_algebra.grothendieck\_basis

<a id="schubmult.rings.free_algebra.grothendieck_basis.GrothendieckBasis"></a>

## GrothendieckBasis Objects

```python
class GrothendieckBasis(FreeAlgebraBasis)
```

Grothendieck basis for FreeAlgebra.

Keys match the Schubert basis shape: ``(Permutation, numvars)``.
The deformation parameter ``beta`` is a class variable so the class can be
passed directly to ``FreeAlgebra`` without requiring a basis constructor.

<a id="schubmult.rings.free_algebra.grothendieck_basis.GrothendieckBasis.product"></a>

#### product

```python
@classmethod
@cache
def product(cls, key1, key2, coeff=S.One)
```

Multiply two keys by transitioning to WordBasis and back.

