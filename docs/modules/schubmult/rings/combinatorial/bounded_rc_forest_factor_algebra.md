<a id="schubmult.rings.combinatorial.bounded_rc_forest_factor_algebra"></a>

# schubmult.rings.combinatorial.bounded\_rc\_forest\_factor\_algebra

<a id="schubmult.rings.combinatorial.bounded_rc_forest_factor_algebra.BoundedRCForestFactorAlgebraElement"></a>

## BoundedRCForestFactorAlgebraElement Objects

```python
class BoundedRCForestFactorAlgebraElement(CrystalGraphRingElement)
```

Element of BoundedRCForestFactorAlgebra: finite linear combinations of Grass tensors.

<a id="schubmult.rings.combinatorial.bounded_rc_forest_factor_algebra.BoundedRCForestFactorAlgebra"></a>

## BoundedRCForestFactorAlgebra Objects

```python
class BoundedRCForestFactorAlgebra(CrystalGraphRing)
```

Tensor-like algebra on tuples of full Grassmannian RC graphs.

Basis keys are tuples (g1, ..., gk) where each gi is a full Grassmannian
RC graph. Simplification rules

<a id="schubmult.rings.combinatorial.bounded_rc_forest_factor_algebra.BoundedRCForestFactorAlgebra.dual_product_on_basis"></a>

#### dual\_product\_on\_basis

```python
def dual_product_on_basis(left_key, right_key)
```

Dual product to the coproduct_on_basis deconcatenation.

<a id="schubmult.rings.combinatorial.bounded_rc_forest_factor_algebra.BoundedRCForestFactorAlgebra.key_to_rc_graph"></a>

#### key\_to\_rc\_graph

```python
def key_to_rc_graph(key) -> RCGraph
```

Evaluate a tensor key to an RCGraph using left-to-right squash_product.

