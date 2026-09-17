<a id="schubmult.rings.combinatorial.grass_tensor_algebra"></a>

# schubmult.rings.combinatorial.grass\_tensor\_algebra

<a id="schubmult.rings.combinatorial.grass_tensor_algebra.GrassTensorAlgebraElement"></a>

## GrassTensorAlgebraElement Objects

```python
class GrassTensorAlgebraElement(CrystalGraphRingElement)
```

Element of GrassTensorAlgebra: finite linear combinations of Grass tensors.

<a id="schubmult.rings.combinatorial.grass_tensor_algebra.GrassTensorAlgebra"></a>

## GrassTensorAlgebra Objects

```python
class GrassTensorAlgebra(CrystalGraphRing)
```

Tensor-like algebra on tuples of full Grassmannian RC graphs.

Basis keys are tuples (g1, ..., gk) where each gi is a full Grassmannian
RC graph. Simplification rules

<a id="schubmult.rings.combinatorial.grass_tensor_algebra.GrassTensorAlgebra.dual_product_on_basis"></a>

#### dual\_product\_on\_basis

```python
def dual_product_on_basis(left_key, right_key)
```

Dual product to the coproduct_on_basis deconcatenation.

<a id="schubmult.rings.combinatorial.grass_tensor_algebra.GrassTensorAlgebra.key_to_rc_graph"></a>

#### key\_to\_rc\_graph

```python
def key_to_rc_graph(key: CrystalGraphTensor | tuple) -> RCGraph
```

Evaluate a tensor key to an RCGraph using left-to-right squash_product.

