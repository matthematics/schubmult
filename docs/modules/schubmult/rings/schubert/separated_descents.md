<a id="schubmult.rings.schubert.separated_descents"></a>

# schubmult.rings.schubert.separated\_descents

Separated-descents ring: Schubert polynomials graded by an explicit number of variables.

`SeparatedDescentsRing` wraps a Schubert-family ring and indexes basis elements by
``(perm, num_vars)`` with ``num_vars >= max_descent(perm)``. The product of
``(u, p)`` and ``(v, q)`` places ``u`` in the first ``p`` variables and ``v`` in the
next ``q``, so descents of the two factors are separated. ``_sep_desc_mul`` computes
this (Samuel) as a twisted ordinary Schubert product by a dominant permutation. Also
carries experimental Pieri/coproduct routines (``pieri_formula``, ``coproduct_test``).

Not to be confused with `schubmult.mult.separated_descents`, which implements the
Fan-Guo-Xiong pipe-puzzle rule for double Grothendieck polynomials.

<a id="schubmult.rings.schubert.separated_descents.complete_sym_positional_perms_down"></a>

#### complete\_sym\_positional\_perms\_down

```python
def complete_sym_positional_perms_down(orig_perm, p, *k, hack_off=None)
```

Descent-side analogue of ``complete_sym_positional_perms``: all ``(perm, degree, sign)`` reachable
from ``orig_perm`` by up to ``p`` Bruhat *descents* swapping a fixed position in ``k`` (1-indexed)
with a still-untouched position. ``hack_off`` bounds the positions considered.

<a id="schubmult.rings.schubert.separated_descents.SeparatedDescentsRing"></a>

## SeparatedDescentsRing Objects

```python
class SeparatedDescentsRing(BaseSchubertRing)
```

Ring with basis ``(perm, num_vars)``; construct with ``SeparatedDescentsRing(Sx.ring)`` and call
as ``ring(perm, num_vars)``. See the module docstring.

<a id="schubmult.rings.schubert.separated_descents.SeparatedDescentsRing.schub_ring"></a>

#### schub\_ring

```python
@property
def schub_ring()
```

The underlying Schubert-family ring used for products.

<a id="schubmult.rings.schubert.separated_descents.SeparatedDescentsRing.pieri_formula"></a>

#### pieri\_formula

```python
def pieri_formula(p, elem)
```

Multiply ``elem`` by the single-row element ``(uncode([p]), 1)`` (adds one variable), via
``complete_sym_positional_perms_down``.

<a id="schubmult.rings.schubert.separated_descents.SeparatedDescentsRing.__init__"></a>

#### \_\_init\_\_

```python
def __init__(ring)
```

Wrap the Schubert-family ``ring`` (inherits its alphabets).

<a id="schubmult.rings.schubert.separated_descents.SeparatedDescentsRing.coproduct_test"></a>

#### coproduct\_test

```python
def coproduct_test(key)
```

Experimental coproduct of the basis element ``key = (perm, num_vars)``, computed by
triangular peeling of the leading code entry via ``pieri_formula``/``_single_coprod_test``.

<a id="schubmult.rings.schubert.separated_descents.SeparatedDescentsRing.mul"></a>

#### mul

```python
def mul(elem1, elem2)
```

Ring product: scalars scale; otherwise each pair of basis elements multiplies via ``_sep_desc_mul``
and lands in degree ``deg1 + deg2`` (terms whose descents exceed that are dropped).

<a id="schubmult.rings.schubert.separated_descents.SeparatedDescentsRing.printing_term"></a>

#### printing\_term

```python
def printing_term(k)
```

The ``SepDescSchubPoly`` display symbol for ``k = (perm, num_vars)``.

<a id="schubmult.rings.schubert.separated_descents.SeparatedDescentsRing.new"></a>

#### new

```python
def new(perm, deg=0)
```

Build ``(perm, deg)`` with ``deg`` raised to at least ``perm``'s last descent; a non-permutation
``perm`` is expanded in the underlying ring and each term given its minimal degree.

<a id="schubmult.rings.schubert.separated_descents.SeparatedDescentsRingElement"></a>

## SeparatedDescentsRingElement Objects

```python
class SeparatedDescentsRingElement(BaseSchubertElement)
```

An element of a `SeparatedDescentsRing`: ``{(perm, num_vars): coeff}``.

<a id="schubmult.rings.schubert.separated_descents.SeparatedDescentsRingElement.coproduct_test"></a>

#### coproduct\_test

```python
def coproduct_test()
```

Experimental coproduct; see `SeparatedDescentsRing.coproduct_test`.

<a id="schubmult.rings.schubert.separated_descents.SeparatedDescentsRingElement.as_ordered_terms"></a>

#### as\_ordered\_terms

```python
def as_ordered_terms(*_, **__)
```

Terms sorted by permutation length, permutation, then ``num_vars`` (sympy printing hook).

