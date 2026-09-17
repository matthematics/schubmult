<a id="schubmult.rings.polynomial_algebra.glide_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.glide\_poly\_basis

<a id="schubmult.rings.polynomial_algebra.glide_poly_basis.glide_monomials"></a>

#### glide\_monomials

```python
def glide_monomials(key)
```

Monomial expansion of the glide polynomial :math:`\mathcal{G}_{key}`.

Returns a dict mapping each exponent tuple ``v`` (a weak composition of the
same length as ``key``) to the integer coefficient of the monomial
:math:`x^v`, i.e. the number of glides of ``key`` with weight ``v``. The
corresponding power of ``beta`` for the weight ``v`` is
``sum(v) - sum(key)`` (the excess), which is constant across all glides of a
given weight, so it need not be stored explicitly.

<a id="schubmult.rings.polynomial_algebra.glide_poly_basis.glide_product"></a>

#### glide\_product

```python
def glide_product(key1, key2)
```

Structure constants for a product of two glide polynomials.

Implements the Littlewood-Richardson rule of O. Pechenik and D. Searles,
"Decompositions of Grothendieck Polynomials" (arXiv:1611.02545), Theorem
4.9, which expands the product of the glide polynomials indexed by the weak
compositions ``key1`` and ``key2`` in the glide basis:

.. math::

    \mathcal{G}_a \, \mathcal{G}_b
        = \sum_c \beta^{|c| - |a| - |b|} \, g_{a,b}^{c} \, \mathcal{G}_c .

Rather than enumerating the genomic shuffle set directly, we compute the
(uniquely determined) coefficients by expanding the product in monomials and
straightening into the glide basis with the leading-term algorithm from the
proof that the glide polynomials form a basis (Theorem 2.6). Because the
excess of a glide of ``v`` equals ``sum(v) - sum(index)``, the power of
``beta`` is recovered from the total degree and only the positive integer
multiplicities :math:`g_{a,b}^{c}` are returned.

Both compositions are padded with trailing zeros to a common length ``n``;
every key ``c`` in the returned dict is a weak composition of length ``n``.

<a id="schubmult.rings.polynomial_algebra.glide_poly_basis.GlidePolyBasis"></a>

## GlidePolyBasis Objects

```python
class GlidePolyBasis(PolynomialBasis)
```

Glide polynomial basis.

Keys are weak compositions. Glide polynomials provide a
basis that refines Grothendieck polynomials and coarsens monomials, with
an efficient combinatorial product rule.

<a id="schubmult.rings.polynomial_algebra.glide_poly_basis.GlidePolyBasis.to_monoms"></a>

#### to\_monoms

```python
def to_monoms(key)
```

Expand a glide key into a dict of monomial exponent tuples.

<a id="schubmult.rings.polynomial_algebra.glide_poly_basis.GlidePolyBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual free algebra basis class (:class:`GlideBasis`).

<a id="schubmult.rings.polynomial_algebra.glide_poly_basis.GlidePolyBasis.expand"></a>

#### expand

```python
def expand(dct)
```

Expand a glide basis dict into a symbolic polynomial expression.

<a id="schubmult.rings.polynomial_algebra.glide_poly_basis.GlidePolyBasis.transition_monomial"></a>

#### transition\_monomial

```python
def transition_monomial(dct)
```

Transition from glide basis to monomial basis.

<a id="schubmult.rings.polynomial_algebra.glide_poly_basis.GlidePolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from glide basis to *other_basis*.

<a id="schubmult.rings.polynomial_algebra.glide_poly_basis.GlidePolyBasis.product"></a>

#### product

```python
@cache
def product(key1, key2, coeff=S.One)
```

Multiply two glide keys using the glide product rule.

