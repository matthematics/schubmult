<a id="schubmult.rings.schubert.double_grothendieck_ring"></a>

# schubmult.rings.schubert.double\_grothendieck\_ring

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckElement"></a>

## DoubleGrothendieckElement Objects

```python
class DoubleGrothendieckElement(BaseSchubertElement)
```

Element of a DoubleGrothendieckRing, stored as {Permutation: coeff}.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing"></a>

## DoubleGrothendieckRing Objects

```python
class DoubleGrothendieckRing(BaseSchubertRing)
```

Ring of double (K-theoretic) Grothendieck polynomials G_w(x, y).

Elements are stored in the G-basis as ``{Permutation: coeff}``. There is
no direct structure-constant formula for the product implemented here:
instead both factors are expanded into the underlying ``DoubleSchubertRing``
(via ``grothendieck_poly_with_ring``), multiplied there, and the product is
converted back to the G-basis with ``to_groth_with_ring``.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.exp_root"></a>

#### exp\_root

```python
@cache
def exp_root(a, b)
```

Polynomial form of the root ``eps_a - eps_b``, i.e. ``(x^root - 1)/beta``.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.exp_weight"></a>

#### exp\_weight

```python
@cache
def exp_weight(index)
```

Polynomial form of the fundamental weight ``eps_index``.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.exp_character"></a>

#### exp\_character

```python
def exp_character(weight)
```

Polynomial form of the multiplicative character ``x^weight``.

Determined by ``exp_root``: ``x^{eps_i} = 1 + beta * y_i``, so that
``(x^{eps_a - eps_b} - 1)/beta`` reproduces ``exp_root(a, b)``.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.chevalley"></a>

#### chevalley

```python
def chevalley(weight, perm, n=None)
```

Lenart--Postnikov K_T-Chevalley formula for ``e^weight * G_perm``.

``weight`` is an integer vector in the ``epsilon`` basis.

