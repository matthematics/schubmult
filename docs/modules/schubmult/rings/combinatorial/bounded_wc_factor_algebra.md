<a id="schubmult.rings.combinatorial.bounded_wc_factor_algebra"></a>

# schubmult.rings.combinatorial.bounded\_wc\_factor\_algebra

`BoundedWCFactorAlgebra`: the `WCGraph` (K-theoretic) analogue of `BoundedRCFactorAlgebra`,
factoring Grothendieck classes into tuples of full Grassmannian WC graphs.

<a id="schubmult.rings.combinatorial.bounded_wc_factor_algebra.BoundedWCFactorPrintingTerm"></a>

## BoundedWCFactorPrintingTerm Objects

```python
class BoundedWCFactorPrintingTerm(PrintingTerm)
```

Display symbol for a `BoundedWCFactorAlgebra` basis key.

<a id="schubmult.rings.combinatorial.bounded_wc_factor_algebra.BoundedWCFactorAlgebraElement"></a>

## BoundedWCFactorAlgebraElement Objects

```python
class BoundedWCFactorAlgebraElement(CrystalGraphRingElement)
```

Element of BoundedWCFactorAlgebra: finite linear combinations of Grass tensors.

<a id="schubmult.rings.combinatorial.bounded_wc_factor_algebra.BoundedWCFactorAlgebra"></a>

## BoundedWCFactorAlgebra Objects

```python
class BoundedWCFactorAlgebra(CrystalGraphRing)
```

Tensor-like algebra on tuples of full Grassmannian WC graphs.

Basis keys are tuples (g1, ..., gk) where each gi is a full Grassmannian
WC graph. Simplification rules

<a id="schubmult.rings.combinatorial.bounded_wc_factor_algebra.BoundedWCFactorAlgebra.key_to_wc_graph"></a>

#### key\_to\_wc\_graph

```python
def key_to_wc_graph(key) -> WCGraph
```

Evaluate a tensor key to an WCGraph using left-to-right squash_product.

