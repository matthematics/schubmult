<a id="schubmult.rings.schubert.double_grothendieck_ring"></a>

# schubmult.rings.schubert.double\_grothendieck\_ring

Double (equivariant K-theoretic) Grothendieck polynomial ring: the ``DGx`` interface.

`DoubleGrothendieckRing` represents ``G_w(x; y)`` with deformation parameter
``beta``. Products go through `schubmult.mult.groth_double.grothmult_double`
(the K-theoretic Monk/Pieri machinery) where available, with a fallback that
expands into the underlying `DoubleSchubertRing`. The ring also exposes the
localization/vanishing data used to convert between the Schubert and
Grothendieck bases (``permuted_subs_dict``, ``product_of_roots``,
``exp_root``, ``chevalley``).

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckElement"></a>

## DoubleGrothendieckElement Objects

```python
class DoubleGrothendieckElement(BaseSchubertElement)
```

Element of a DoubleGrothendieckRing, stored as {Permutation: coeff}.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckElement.as_polynomial"></a>

#### as\_polynomial

```python
def as_polynomial()
```

Expand to an explicit polynomial: ``sum coeff * G_w(x; y)``.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckElement.perm_subs"></a>

#### perm\_subs

```python
def perm_subs(perm)
```

Localize at the torus fixed point ``perm``: substitute ``x_i -> (-) y_{perm(i)}`` (formal inverse).

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckElement.simplify"></a>

#### simplify

```python
def simplify(factor=True)
```

Return a copy with each coefficient put in cancelled (and, by default, factored) rational
normal form in ``y`` and ``beta``, dropping terms whose coefficient simplifies to zero.

Products in this ring leave coefficients as unsimplified rational expressions; this
makes them readable, e.g. ``(y_1 - y_2)/(1 + beta*y_2)``.

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

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.perm_subs"></a>

#### perm\_subs

```python
def perm_subs(elem, perm)
```

Localize ``elem`` at ``perm``: expand into double Schubert polynomials and substitute
``x_i -> -y_{perm(i)} / (1 + beta y_{perm(i)})``.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.double_mul"></a>

#### double\_mul

```python
@property
def double_mul()
```

`schubmult.mult.groth_double.grothmult_double`.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.single_mul"></a>

#### single\_mul

```python
@property
def single_mul()
```

`schubmult.mult.groth.grothmult_py`.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.beta"></a>

#### beta

```python
@property
def beta()
```

The deformation parameter.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.single_variable"></a>

#### single\_variable

```python
@property
def single_variable()
```

`schubmult.mult.groth_double.single_variable_groth` (K-theoretic Monk rule for ``x_k``).

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.vanish_subs_dict"></a>

#### vanish\_subs\_dict

```python
@cached_property
def vanish_subs_dict()
```

Localization at the identity: ``x_i -> -y_i / (1 + beta y_i)`` for the first 50 variables.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.permuted_subs_dict"></a>

#### permuted\_subs\_dict

```python
def permuted_subs_dict(perm, length=None)
```

Localization at ``perm``: ``x_i -> -y_{perm(i)} / (1 + beta y_{perm(i)})`` for ``i <= length``.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.product_of_roots"></a>

#### product\_of\_roots

```python
@cache
def product_of_roots(perm)
```

``prod_{(a,b) in Inv(perm^-1)} (y_a (-) y_b)``: the localization of ``G_perm`` at itself (Euler class).

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

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.div_by_product_of_roots"></a>

#### div\_by\_product\_of\_roots

```python
def div_by_product_of_roots(expr, perm)
```

Divide ``expr`` by ``product_of_roots(perm)`` one linear factor at a time, leaving any
non-exact factors in the denominator (so the result may be a rational function).

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.mult_poly_double"></a>

#### mult\_poly\_double

```python
@property
def mult_poly_double()
```

`schubmult.mult.groth_double.mult_poly_groth_double`.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.mult_poly_single"></a>

#### mult\_poly\_single

```python
@property
def mult_poly_single()
```

`schubmult.mult.groth.mult_poly_groth`.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.schub_as_groth"></a>

#### schub\_as\_groth

```python
@cache
def schub_as_groth(perm)
```

The double Schubert polynomial ``S_perm`` expanded in the Grothendieck basis (cached).

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.from_double_schubert_elem"></a>

#### from\_double\_schubert\_elem

```python
def from_double_schubert_elem(elem)
```

Convert a `DoubleSchubertElement` into this ring's Grothendieck basis.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.mul_expr"></a>

#### mul\_expr

```python
def mul_expr(elem, x)
```

Multiply by an expression: single ``x`` variables via the K-Monk rule, ``Add``/``Mul``/``Pow``
recursively, anything else as a coefficient.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.mul"></a>

#### mul

```python
def mul(elem, other)
```

Ring product via ``_best_effort_grothmult_double``.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.from_expr"></a>

#### from\_expr

```python
def from_expr(expr)
```

Convert a polynomial into the Grothendieck basis by multiplying the identity by it.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.cached_schubpoly"></a>

#### cached\_schubpoly

```python
@cache
def cached_schubpoly(k)
```

The explicit ``G_k(x; y)``, as a sum of `WCGraph` monomials weighted by ``beta^excess``.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.printing_term"></a>

#### printing\_term

```python
def printing_term(k, prefix="")
```

The ``DoubleGrothendieckPoly`` display symbol for basis element ``k``.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.new"></a>

#### new

```python
def new(x)
```

Build an element from a permutation/Lehmer list, an element of this ring, or a polynomial expression.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.from_dict"></a>

#### from\_dict

```python
def from_dict(dct)
```

Build an element from ``{Permutation: coeff}``, dropping exact zeros.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DGx"></a>

#### DGx

```python
def DGx(x, genset=GeneratingSet("y"))
```

Construct a double Grothendieck element in ``x`` with coefficient alphabet ``genset`` (a
`GeneratingSet`, a label string, or ``"0"`` for the zero alphabet).

