<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring"></a>

# schubmult.rings.schubert.parabolic\_quantum\_double\_schubert\_ring

Parabolic quantum double Schubert polynomial ring: the ``QPDSx`` interface.

`ParabolicQuantumDoubleSchubertRing` models the quantum cohomology of a partial
flag variety with block sizes ``index_comp``. Basis permutations must be
parabolic (increasing within each block). Products are computed in the full
flag quantum ring and projected down via the Peterson-Woodward comparison
(`schubmult.mult.quantum_double.apply_peterson_woodward`, through
``process_coeff_dict``). The parabolic quantum elementary symmetric functions
acquire a ``q``-correction at each block boundary.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertElement"></a>

## ParabolicQuantumDoubleSchubertElement Objects

```python
class ParabolicQuantumDoubleSchubertElement(BaseSchubertElement)
```

An element of a `ParabolicQuantumDoubleSchubertRing`.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertElement.index_comp"></a>

#### index\_comp

```python
@property
def index_comp()
```

The ring's block-size composition.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertElement.kill_ideal"></a>

#### kill\_ideal

```python
def kill_ideal()
```

Drop basis permutations longer than ``sum(index_comp)`` (those lie in the ideal cut out by the parabolic).

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing"></a>

## ParabolicQuantumDoubleSchubertRing Objects

```python
class ParabolicQuantumDoubleSchubertRing(BaseSchubertRing)
```

Quantum double Schubert polynomials for the partial flag variety with block sizes ``index_comp``.
Construct via ``QPDSx(*index_comp)([perm])`` or ``make_parabolic_quantum_basis``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.__init__"></a>

#### \_\_init\_\_

```python
def __init__(genset, coeff_genset, index_comp)
```

**Arguments**:

- `genset` - Primary ``x`` alphabet.
- `coeff_genset` - Coefficient ``y`` alphabet.
- `index_comp` - Composition of block sizes; ``sum(index_comp)`` is the ambient ``n``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.symbol_elem_func"></a>

#### symbol\_elem\_func

```python
@property
def symbol_elem_func()
```

Symbolic elementary symmetric function ``e_p_k`` combined with complete symmetric corrections in ``-y``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.elem_sym_subs"></a>

#### elem\_sym\_subs

```python
def elem_sym_subs(kk)
```

Substitution dict ``{e_p_k: elem_sym(p, k, x, 0)}`` for all ``1 <= p <= k <= kk``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.parabolic_index"></a>

#### parabolic\_index

```python
@property
def parabolic_index()
```

1-indexed positions of the simple reflections inside the parabolic subgroup (within-block positions).

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.quantum_basis"></a>

#### quantum\_basis

```python
@property
def quantum_basis()
```

The full-flag `QuantumDoubleSchubertRing` over the same alphabets.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.classical_basis"></a>

#### classical\_basis

```python
@property
def classical_basis()
```

The classical `DoubleSchubertRing` over the same alphabets.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.elem_sym"></a>

#### elem\_sym

```python
def elem_sym(p, k, varl1, varl2)
```

Parabolic quantum double elementary symmetric polynomial ``E_p(x_1..x_k; y)``: classical below
the first block boundary, with a ``q_j``-correction at each block boundary ``N_j``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.index_comp"></a>

#### index\_comp

```python
@property
def index_comp()
```

The block-size composition.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.process_coeff_dict"></a>

#### process\_coeff\_dict

```python
def process_coeff_dict(coeff_dict)
```

Project a full-flag quantum coefficient dict onto this parabolic ring via Peterson-Woodward,
extending the parabolic index if any permutation exceeds the ambient ``n``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.cached_product"></a>

#### cached\_product

```python
@cache
def cached_product(u, v, basis2)
```

Full-flag quantum double product (generic alphabets, substituted) then projected via ``process_coeff_dict``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.in_quantum_basis"></a>

#### in\_quantum\_basis

```python
def in_quantum_basis(elem)
```

Expand into the full-flag quantum double Schubert basis via ``quantum_elem_func``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.in_classical_basis"></a>

#### in\_classical\_basis

```python
def in_classical_basis(elem)
```

Expand into the classical double Schubert basis via ``quantum_as_classical_schubpoly``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.classical_in_basis"></a>

#### classical\_in\_basis

```python
@cache
def classical_in_basis(k)
```

Express the classical ``S_k`` in this parabolic quantum basis, by iteratively subtracting
off lower-order corrections until the polynomials agree.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.classical_elem_func"></a>

#### classical\_elem\_func

```python
@property
def classical_elem_func()
```

Parabolic quantum elementary symmetric function valued in the classical `DoubleSchubertRing`.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.quantum_elem_func"></a>

#### quantum\_elem\_func

```python
@property
def quantum_elem_func()
```

Parabolic quantum elementary symmetric function valued in the full-flag `QuantumDoubleSchubertRing`.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.printing_term"></a>

#### printing\_term

```python
def printing_term(k)
```

The ``PQDSchubPoly`` display symbol for basis element ``k``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.quantum_as_classical_schubpoly"></a>

#### quantum\_as\_classical\_schubpoly

```python
@cache
def quantum_as_classical_schubpoly(perm)
```

``S^{q,P}_perm`` expanded in the classical double Schubert basis, against the appropriate longest element.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.cached_schubpoly"></a>

#### cached\_schubpoly

```python
@cache
def cached_schubpoly(k)
```

The explicit parabolic quantum double Schubert polynomial for ``k``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.cached_positive_product"></a>

#### cached\_positive\_product

```python
@cache
def cached_positive_product(u, v, basis2)
```

Positive variant of ``cached_product`` via ``schubmult_q_generic_partial_posify``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.double_mul"></a>

#### double\_mul

```python
@property
def double_mul()
```

``schubmult_q_double_fast`` followed by ``process_coeff_dict``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.single_mul"></a>

#### single\_mul

```python
@property
def single_mul()
```

``schubmult_q_fast`` followed by ``process_coeff_dict``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.mult_poly_single"></a>

#### mult\_poly\_single

```python
@property
def mult_poly_single()
```

`schubmult.mult.quantum.mult_poly_q`.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.mult_poly_double"></a>

#### mult\_poly\_double

```python
@property
def mult_poly_double()
```

`schubmult.mult.quantum_double.mult_poly_q_double`.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.from_expr"></a>

#### from\_expr

```python
def from_expr(expr)
```

Convert a polynomial to this basis by peeling off leading monomials; raises ``ValueError`` if
``expr`` lacks the within-block symmetry the parabolic ring requires.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.mul_expr"></a>

#### mul\_expr

```python
def mul_expr(elem, x)
```

Multiply by an expression by first converting it into this basis.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.__call__"></a>

#### \_\_call\_\_

```python
def __call__(x)
```

Build an element from a parabolic permutation/Lehmer list or an expression; raises ``ValueError``
if the permutation is not parabolic for this ring's blocks.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.make_parabolic_quantum_basis"></a>

#### make\_parabolic\_quantum\_basis

```python
def make_parabolic_quantum_basis(index_comp, coeff_genset)
```

The `ParabolicQuantumDoubleSchubertRing` in ``x`` for block sizes ``index_comp`` and coefficient alphabet ``coeff_genset``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.QPDSx_index"></a>

#### QPDSx\_index

```python
def QPDSx_index(*args)
```

Return a constructor ``f(x, coeff_genset="y")`` building parabolic quantum double elements for block sizes ``args``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.QPDSx"></a>

#### QPDSx

```python
@cache
def QPDSx(*args)
```

Cached constructor for block sizes ``args``; e.g. ``QPDSx(2, 1)([2, 1, 3])``.

