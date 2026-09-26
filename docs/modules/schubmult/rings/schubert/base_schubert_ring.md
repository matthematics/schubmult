<a id="schubmult.rings.schubert.base_schubert_ring"></a>

# schubmult.rings.schubert.base\_schubert\_ring

Abstract base classes shared by every Schubert-family ring.

`BaseSchubertRing` holds a primary generating set (``genset``, the ``x`` variables)
and a coefficient generating set (``coeff_genset``, the ``y``/``z`` variables, or
``None`` for single Schubert polynomials), and defines the ring-level hooks that
concrete rings fill in: which multiplication kernel to use, how to expand a basis
element to a polynomial, how to print it, and how to change basis. `BaseSchubertElement`
is the corresponding dict-like element type (``{Permutation: coefficient}``).

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertElement"></a>

## BaseSchubertElement Objects

```python
class BaseSchubertElement(BaseRingElement)
```

A linear combination of Schubert-family basis elements, stored as ``{Permutation: coeff}``.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertElement.mult_poly"></a>

#### mult\_poly

```python
def mult_poly(poly)
```

Multiply this element by an arbitrary polynomial ``poly`` in the ring's ``genset`` variables.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertElement.in_schubert_schur_basis"></a>

#### in\_schubert\_schur\_basis

```python
def in_schubert_schur_basis(numvars)
```

Expand into the Schubert-tensor-Schur basis of the tensor square ring, splitting off the
symmetric part in the last ``numvars`` variables.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertElement.in_SEM_basis"></a>

#### in\_SEM\_basis

```python
def in_SEM_basis(elem_func=None)
```

Expand as a polynomial in elementary symmetric functions (the "SEM" presentation), using
``elem_func`` (default: the ring's symbolic ``symbol_elem_func``) as the elementary symmetric symbol.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertElement._repr_latex_"></a>

#### \_repr\_latex\_

```python
def _repr_latex_()
```

Disabled so notebooks show the fast text form; use ``latex(elem)`` or ``pretty(elem)`` explicitly.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertElement.as_ordered_terms"></a>

#### as\_ordered\_terms

```python
def as_ordered_terms(*_, **__)
```

Terms ``coeff * basis_symbol`` sorted by permutation length then lexicographically (sympy printing hook).

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertElement.expand"></a>

#### expand

```python
def expand(deep=True, *args, **kwargs)
```

With ``deep=True`` (default) expand to an explicit polynomial in the variables; with
``deep=False`` only expand each coefficient, keeping the Schubert basis.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertElement.as_expr"></a>

#### as\_expr

```python
def as_expr()
```

Sum of the ``as_terms()`` as a sympy ``Add``.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertElement.as_polynomial"></a>

#### as\_polynomial

```python
def as_polynomial()
```

Expand to an explicit polynomial: ``sum coeff * SchubertPoly(perm)``.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertElement.as_classical"></a>

#### as\_classical

```python
def as_classical()
```

Re-express in the classical (non-quantum) Schubert basis.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertElement.as_quantum"></a>

#### as\_quantum

```python
def as_quantum()
```

Re-express in the quantum Schubert basis.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertElement.almosteq"></a>

#### almosteq

```python
def almosteq(other)
```

Equality up to coefficient expansion (handles elements of different but compatible rings).

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertElement.strip_zeros"></a>

#### strip\_zeros

```python
def strip_zeros()
```

Drop basis elements whose coefficient is exactly zero.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing"></a>

## BaseSchubertRing Objects

```python
class BaseSchubertRing(BaseRing)
```

Abstract base ring for Schubert-family polynomials.

Concrete subclasses supply the multiplication kernels (``double_mul``/``single_mul``,
``mult_poly_single``/``mult_poly_double``), the basis-element expansion
(``cached_schubpoly``), printing (``printing_term``), coercion, and basis changes.
Two rings compare equal iff they have the same type and generating sets.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.__init__"></a>

#### \_\_init\_\_

```python
def __init__(genset, coeff_genset, domain=None)
```

**Arguments**:

- `genset` - Primary generating set (the ``x`` variables).
- `coeff_genset` - Coefficient generating set (``y``/``z``), or a set with ``label=None`` for single rings.
- `domain` - Optional coefficient domain passed to `BaseRing`.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.mul"></a>

#### mul

```python
def mul(elem, other)
```

Multiply two elements via `_mul_schub_dicts`, which dispatches to the appropriate
`schubmult.mult` kernel based on both rings' generating sets; scalars go through `BaseRing.mul`.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.new"></a>

#### new

```python
def new(x)
```

Hook: build an element from ``x`` (permutation, code, or expression).

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.printing_term"></a>

#### printing\_term

```python
def printing_term(k)
```

Hook: the sympy symbol displayed for basis element ``k``.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.coproduct_on_basis"></a>

#### coproduct\_on\_basis

```python
def coproduct_on_basis(k)
```

Hook: coproduct of basis element ``k`` in the tensor square ring.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.is_elem_mul_type"></a>

#### is\_elem\_mul\_type

```python
def is_elem_mul_type(elem)
```

Hook: whether ``elem`` should be multiplied via the elementary-symmetric fast path.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.elem_mul"></a>

#### elem\_mul

```python
def elem_mul(ring_elem, elem)
```

Hook: elementary-symmetric fast-path multiplication.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.elem_sym"></a>

#### elem\_sym

```python
@property
def elem_sym()
```

Hook: the elementary symmetric polynomial function used by this ring.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.symbol_elem_func"></a>

#### symbol\_elem\_func

```python
@property
def symbol_elem_func()
```

Hook: symbolic (unevaluated) elementary symmetric function for ``in_SEM_basis``.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.elem_sym_subs"></a>

#### elem\_sym\_subs

```python
def elem_sym_subs(kk)
```

Hook: substitution dict turning the symbolic elementary symmetric symbols back into polynomials.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.domain_new"></a>

#### domain\_new

```python
def domain_new(element, orig_domain=None)
```

Coerce ``element`` into the coefficient domain, refusing anything containing a ``genset`` variable.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.genset"></a>

#### genset

```python
@property
def genset()
```

Primary generating set (the ``x`` variables).

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.coeff_genset"></a>

#### coeff\_genset

```python
@property
def coeff_genset()
```

Coefficient generating set (``y``/``z`` variables); ``label`` is ``None`` for single rings.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.in_quantum_basis"></a>

#### in\_quantum\_basis

```python
def in_quantum_basis(elem)
```

Hook: re-express ``elem`` in the quantum Schubert basis.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.in_classical_basis"></a>

#### in\_classical\_basis

```python
def in_classical_basis(elem)
```

Hook: re-express ``elem`` in the classical Schubert basis.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.quantum_schubpoly"></a>

#### quantum\_schubpoly

```python
def quantum_schubpoly(perm)
```

Hook: the quantum Schubert polynomial for ``perm``.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.cached_product"></a>

#### cached\_product

```python
def cached_product(u, v, basis2)
```

Hook: cached structure constants of ``S_u * S_v`` with ``v`` in ``basis2``.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.cached_positive_product"></a>

#### cached\_positive\_product

```python
def cached_positive_product(u, v, basis2)
```

Hook: like ``cached_product`` but with manifestly positive coefficients.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.mul_expr"></a>

#### mul\_expr

```python
def mul_expr(elem, x)
```

Hook: multiply ``elem`` by a symbolic expression ``x``.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.double_mul"></a>

#### double\_mul

```python
@property
def double_mul()
```

Hook: the double-Schubert multiplication kernel (e.g. ``schubmult_double``).

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.single_mul"></a>

#### single\_mul

```python
@property
def single_mul()
```

Hook: the single-Schubert multiplication kernel (e.g. ``schubmult_py``).

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.mult_poly_single"></a>

#### mult\_poly\_single

```python
@property
def mult_poly_single()
```

Hook: the single-variant multiply-by-polynomial kernel.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.mult_poly_double"></a>

#### mult\_poly\_double

```python
@property
def mult_poly_double()
```

Hook: the double-variant multiply-by-polynomial kernel.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.quantum_elem_func"></a>

#### quantum\_elem\_func

```python
@property
def quantum_elem_func()
```

Hook: the quantum elementary symmetric function.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.cached_schubpoly"></a>

#### cached\_schubpoly

```python
def cached_schubpoly(k)
```

Hook: the (cached) explicit polynomial for basis element ``k``.

