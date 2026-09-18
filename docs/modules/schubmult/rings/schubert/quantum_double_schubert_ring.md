<a id="schubmult.rings.schubert.quantum_double_schubert_ring"></a>

# schubmult.rings.schubert.quantum\_double\_schubert\_ring

Quantum double Schubert polynomial ring: the ``QDSx`` interface.

`QuantumDoubleSchubertRing` represents quantum double Schubert polynomials
``S^q_w(x; y)`` with quantum parameters ``q_1, q_2, ...`` (the module-level
``q_var``), dispatching products to `schubmult.mult.quantum_double`. Elements
can be converted to/from the classical basis via ``as_classical``/``as_quantum``.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.is_fact_elem_sym"></a>

#### is\_fact\_elem\_sym

```python
def is_fact_elem_sym(obj)
```

Whether ``obj`` is an (unevaluated) quantum factorial elementary symmetric function.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertElement"></a>

## QuantumDoubleSchubertElement Objects

```python
class QuantumDoubleSchubertElement(BaseSchubertElement)
```

An element of a `QuantumDoubleSchubertRing`: ``{Permutation: coeff}`` in the basis ``S^q_w(x; y)``.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertElement.subs"></a>

#### subs

```python
def subs(old, new)
```

Substitute by round-tripping through the classical basis.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing"></a>

## QuantumDoubleSchubertRing Objects

```python
class QuantumDoubleSchubertRing(BaseSchubertRing)
```

The ring of quantum double Schubert polynomials; ``QDSx`` is the standard constructor.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.printing_term"></a>

#### printing\_term

```python
def printing_term(k)
```

The ``QDSchubPoly`` display symbol for basis element ``k``.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.symbol_elem_func"></a>

#### symbol\_elem\_func

```python
@property
def symbol_elem_func()
```

`QFactorialElemSym`.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.elem_sym_subs"></a>

#### elem\_sym\_subs

```python
def elem_sym_subs(kk)
```

Substitution dict ``{e_p_k: elem_sym_poly_q(p, k, x)}`` for all ``1 <= p <= k <= kk``.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.elem_sym"></a>

#### elem\_sym

```python
@property
def elem_sym()
```

`QFactorialElemSym`.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.is_elem_mul_type"></a>

#### is\_elem\_mul\_type

```python
def is_elem_mul_type(other)
```

Whether ``other`` is a quantum factorial elementary symmetric function.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.elem_mul"></a>

#### elem\_mul

```python
def elem_mul(ring_elem, elem)
```

Multiply by a quantum factorial elementary symmetric function via the quantum positional
Pieri rule (``elem_sym_positional_perms_q``), each term carrying its ``q``-monomial.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.cached_product"></a>

#### cached\_product

```python
@cache
def cached_product(u, v, basis2)
```

Structure constants via ``schubmult_q_double_fast``.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.in_quantum_basis"></a>

#### in\_quantum\_basis

```python
def in_quantum_basis(elem)
```

Identity (this ring is already quantum).

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.in_classical_basis"></a>

#### in\_classical\_basis

```python
def in_classical_basis(elem)
```

Expand each ``S^q_w`` as a classical double Schubert element via ``quantum_as_classical_schubpoly``.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.classical_elem_func"></a>

#### classical\_elem\_func

```python
@property
def classical_elem_func()
```

Quantum elementary symmetric function valued in the classical `DoubleSchubertRing`, via the
recursion ``E_p(k) = (x_k - y_{k-p+1}) E_{p-1}(k-1) + E_p(k-1) + q_{k-1} E_{p-2}(k-2)``.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.quantum_as_classical_schubpoly"></a>

#### quantum\_as\_classical\_schubpoly

```python
@cache
def quantum_as_classical_schubpoly(perm)
```

``S^q_perm`` expanded in the classical double Schubert basis (cached).

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.cached_schubpoly"></a>

#### cached\_schubpoly

```python
@cache
def cached_schubpoly(k)
```

The explicit quantum double Schubert polynomial ``S^q_k(x; y)`` (cached).

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.cached_positive_product"></a>

#### cached\_positive\_product

```python
@cache
def cached_positive_product(u, v, basis2)
```

Structure constants with (partially) positive coefficients via ``schubmult_q_generic_partial_posify``.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.double_mul"></a>

#### double\_mul

```python
@property
def double_mul()
```

`schubmult.mult.quantum_double.schubmult_q_double_fast`.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.single_mul"></a>

#### single\_mul

```python
@property
def single_mul()
```

`schubmult.mult.quantum.schubmult_q_fast`.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.mult_poly_single"></a>

#### mult\_poly\_single

```python
@property
def mult_poly_single()
```

`schubmult.mult.quantum.mult_poly_q`.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.positive_elem_sym_rep"></a>

#### positive\_elem\_sym\_rep

```python
def positive_elem_sym_rep(perm, index=1)
```

Manifestly positive expansion of ``S^q_perm`` in quantum factorial elementary symmetric functions (forward).

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.positive_elem_sym_rep_backward"></a>

#### positive\_elem\_sym\_rep\_backward

```python
def positive_elem_sym_rep_backward(perm)
```

Like ``positive_elem_sym_rep`` but peeling from the last descent backward.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.mult_poly_double"></a>

#### mult\_poly\_double

```python
@property
def mult_poly_double()
```

`schubmult.mult.quantum_double.mult_poly_q_double`.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.from_expr"></a>

#### from\_expr

```python
def from_expr(expr)
```

Convert a polynomial to the quantum basis by repeatedly peeling off the leading monomial's
Schubert term; falls back to ``mul_expr`` on the identity if that fails.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.handle_sympoly"></a>

#### handle\_sympoly

```python
def handle_sympoly(other)
```

Evaluate symmetric-function coefficients to polynomials.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.mul_expr"></a>

#### mul\_expr

```python
def mul_expr(elem, x)
```

Multiply by an expression: single ``x`` variables via ``mult_poly_q_double``, quantum factorial
elementary symmetric functions via ``elem_mul``, ``Add``/``Mul``/``Pow`` recursively, else as a coefficient.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.new"></a>

#### new

```python
def new(x)
```

Build an element from a permutation/Lehmer list or a polynomial expression.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QDSx"></a>

#### QDSx

```python
def QDSx(x, genset=GeneratingSet("y"))
```

Construct a quantum double Schubert element in ``x`` with coefficient alphabet ``genset``
(a `GeneratingSet` or a label string); e.g. ``QDSx([3, 1, 2])`` is ``S^q_{312}(x; y)``.

