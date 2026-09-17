<a id="schubmult.rings.schubert.grothendieck_ring"></a>

# schubmult.rings.schubert.grothendieck\_ring

<a id="schubmult.rings.schubert.grothendieck_ring.GrothendieckElement"></a>

## GrothendieckElement Objects

```python
class GrothendieckElement(BaseSchubertElement)
```

Element of a GrothendieckRing, stored as {Permutation: coeff}.

<a id="schubmult.rings.schubert.grothendieck_ring.GrothendieckRing"></a>

## GrothendieckRing Objects

```python
class GrothendieckRing(BaseSchubertRing)
```

Ring of Grothendieck polynomials.

A deformation of the Schubert polynomial ring with parameter beta.
Basis elements G_w satisfy G_u * G_v = sum_w c^w_{u,v}(beta) G_w
where c^w_{u,v}(beta) are polynomials in beta with integer coefficients.
When beta=0, recovers ordinary Schubert polynomials.

Parameters
----------
genset : GeneratingSet
    The generating set (variable alphabet).
beta : sympy/symengine symbol, optional
    The deformation parameter. Defaults to Symbol("β").

