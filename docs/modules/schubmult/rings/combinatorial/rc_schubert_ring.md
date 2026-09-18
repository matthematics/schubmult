<a id="schubmult.rings.combinatorial.rc_schubert_ring"></a>

# schubmult.rings.combinatorial.rc\_schubert\_ring

Orphaned earlier variant of `bounded_rc_factor_algebra` (defines the same class names, is not
imported anywhere, and is not exported from the package). Kept for reference; prefer
`schubmult.rings.combinatorial.bounded_rc_factor_algebra`.

<a id="schubmult.rings.combinatorial.rc_schubert_ring.BoundedRCFactorPrintingTerm"></a>

## BoundedRCFactorPrintingTerm Objects

```python
class BoundedRCFactorPrintingTerm(PrintingTerm)
```

Display symbol for a basis key (orphaned variant).

<a id="schubmult.rings.combinatorial.rc_schubert_ring.BoundedRCFactorAlgebraElement"></a>

## BoundedRCFactorAlgebraElement Objects

```python
class BoundedRCFactorAlgebraElement(CrystalGraphRingElement)
```

Element of BoundedRCFactorAlgebra: finite linear combinations of Grass tensors.

<a id="schubmult.rings.combinatorial.rc_schubert_ring.BoundedRCFactorAlgebraElement.prune"></a>

#### prune

```python
def prune()
```

Merge terms whose evaluated RC graphs share forest/omega invariants.

Terms are grouped by
``(rc.forest_weight, rc.omega_invariant[1])`` where
``rc = self.ring.key_to_rc_graph(key)``.
One key per group is kept (first encountered), and coefficients are summed.

<a id="schubmult.rings.combinatorial.rc_schubert_ring.BoundedRCFactorAlgebra"></a>

## BoundedRCFactorAlgebra Objects

```python
class BoundedRCFactorAlgebra(CrystalGraphRing)
```

Tensor-like algebra on tuples of full Grassmannian RC graphs.

Basis keys are tuples (g1, ..., gk) where each gi is a full Grassmannian
RC graph. Simplification rules

<a id="schubmult.rings.combinatorial.rc_schubert_ring.BoundedRCFactorAlgebra.dual_product_on_basis"></a>

#### dual\_product\_on\_basis

```python
def dual_product_on_basis(left_key, right_key)
```

Dual product to the coproduct_on_basis deconcatenation.

<a id="schubmult.rings.combinatorial.rc_schubert_ring.BoundedRCFactorAlgebra.key_to_rc_graph"></a>

#### key\_to\_rc\_graph

```python
def key_to_rc_graph(key) -> RCGraph
```

Evaluate a tensor key to an RCGraph using left-to-right squash_product.

