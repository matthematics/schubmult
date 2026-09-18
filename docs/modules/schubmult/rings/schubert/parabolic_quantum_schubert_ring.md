<a id="schubmult.rings.schubert.parabolic_quantum_schubert_ring"></a>

# schubmult.rings.schubert.parabolic\_quantum\_schubert\_ring

Parabolic quantum (single) Schubert polynomial ring: the ``QPSx`` interface.

`ParabolicQuantumSingleSchubertRing` is a `ParabolicQuantumDoubleSchubertRing`
with a zero coefficient alphabet, indexed by a composition ``index_comp``
specifying the parabolic subgroup's block sizes.

<a id="schubmult.rings.schubert.parabolic_quantum_schubert_ring.ParabolicQuantumSingleSchubertRing"></a>

## ParabolicQuantumSingleSchubertRing Objects

```python
class ParabolicQuantumSingleSchubertRing(ParabolicQuantumDoubleSchubertRing)
```

Quantum Schubert polynomials for the partial flag variety with block sizes ``index_comp``.
Construct via ``QPSx(*index_comp)``; basis permutations must be parabolic for the given blocks.

<a id="schubmult.rings.schubert.parabolic_quantum_schubert_ring.ParabolicQuantumSingleSchubertRing.cached_schubpoly"></a>

#### cached\_schubpoly

```python
@cache
def cached_schubpoly(k)
```

The explicit parabolic quantum Schubert polynomial for ``k``, expanded against the appropriate
longest element (extended if ``k`` is larger than the ring's default).

<a id="schubmult.rings.schubert.parabolic_quantum_schubert_ring.ParabolicQuantumSingleSchubertRing.elem_sym"></a>

#### elem\_sym

```python
def elem_sym(p, k, varl1, varl2)
```

Parabolic quantum elementary symmetric polynomial ``E_p(x_1..x_k)``: classical below the first
block boundary, with a ``q``-correction term at each block boundary ``N_j``.

<a id="schubmult.rings.schubert.parabolic_quantum_schubert_ring.ParabolicQuantumSingleSchubertRing.coeff_genset"></a>

#### coeff\_genset

```python
@property
def coeff_genset()
```

Always the zero alphabet.

<a id="schubmult.rings.schubert.parabolic_quantum_schubert_ring.ParabolicQuantumSingleSchubertRing.cached_product"></a>

#### cached\_product

```python
@cache
def cached_product(u, v, basis2)
```

Full-flag quantum product (``schubmult_q_fast`` or the generic double kernel), then projected
to the parabolic ring via ``process_coeff_dict`` (Peterson-Woodward).

<a id="schubmult.rings.schubert.parabolic_quantum_schubert_ring.ParabolicQuantumSingleSchubertRing.cached_positive_product"></a>

#### cached\_positive\_product

```python
@cache
def cached_positive_product(u, v, basis2)
```

Same as ``cached_product``.

<a id="schubmult.rings.schubert.parabolic_quantum_schubert_ring.ParabolicQuantumSingleSchubertRing.__call__"></a>

#### \_\_call\_\_

```python
def __call__(x)
```

Build an element from a parabolic permutation/Lehmer list or an expression; raises ``ValueError``
if the permutation is not parabolic for this ring's blocks.

<a id="schubmult.rings.schubert.parabolic_quantum_schubert_ring.make_single_parabolic_quantum_basis"></a>

#### make\_single\_parabolic\_quantum\_basis

```python
def make_single_parabolic_quantum_basis(index_comp)
```

The `ParabolicQuantumSingleSchubertRing` in ``x`` for block sizes ``index_comp``.

<a id="schubmult.rings.schubert.parabolic_quantum_schubert_ring.QPSx"></a>

#### QPSx

```python
@cache
def QPSx(*args)
```

Cached `ParabolicQuantumSingleSchubertRing` for block sizes ``args``, e.g. ``QPSx(2, 1)([2, 1, 3])``.

