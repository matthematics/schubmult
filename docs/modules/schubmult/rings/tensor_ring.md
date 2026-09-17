<a id="schubmult.rings.tensor_ring"></a>

# schubmult.rings.tensor\_ring

<a id="schubmult.rings.tensor_ring.TensorRing"></a>

## TensorRing Objects

```python
class TensorRing(BaseRing)
```

<a id="schubmult.rings.tensor_ring.TensorRing.coproduct_on_basis"></a>

#### coproduct\_on\_basis

```python
def coproduct_on_basis(k)
```

Compute coproduct of a basis element in the tensor ring.

For k = (k_1, k_2, ..., k_n) in ring_1 ⊗ ring_2 ⊗ ... ⊗ ring_n,
Δ(k_1 ⊗ k_2 ⊗ ... ⊗ k_n) = (⊗ Δ(k_i))

This properly interlaces the individual coproducts.

<a id="schubmult.rings.tensor_ring.TensorRing.mul"></a>

#### mul

```python
def mul(elem1, elem2)
```

Multiply two elements in the tensor ring.

(a1 ⊗ a2 ⊗ ... ⊗ an) * (b1 ⊗ b2 ⊗ ... ⊗ bn) = (a1*b1) ⊗ (a2*b2) ⊗ ... ⊗ (an*bn)

<a id="schubmult.rings.tensor_ring.TensorRingElement"></a>

## TensorRingElement Objects

```python
class TensorRingElement(BaseRingElement)
```

<a id="schubmult.rings.tensor_ring.TensorRingElement.coproduct"></a>

#### coproduct

```python
def coproduct()
```

Override coproduct to use the correct target ring.

