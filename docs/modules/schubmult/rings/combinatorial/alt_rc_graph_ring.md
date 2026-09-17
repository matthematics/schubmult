<a id="schubmult.rings.combinatorial.alt_rc_graph_ring"></a>

# schubmult.rings.combinatorial.alt\_rc\_graph\_ring

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement"></a>

## AltRCGraphRingElement Objects

```python
class AltRCGraphRingElement(CrystalGraphRingElement,
                            SchubertMonomialRingElement)
```

AltRCGraphRing elements are linear combinations of RCGraph basis elements.

The product % is the polynomial product. Currently only defined when the right side
is a dominant RC graph.

The Leibniz rule should hold for % somehow. Claude's idea is to define the ambiguous term in the Leibniz formula instead of trying
to do this directly.

The product * is well defined for any pair of RC graphs and is the dual product.

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.__mod__"></a>

#### \_\_mod\_\_

```python
def __mod__(other)
```

Polynomial product: self % other.
Currently only defined when `other` is a dominant RC graph.

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.divdiff_perm"></a>

#### divdiff\_perm

```python
def divdiff_perm(perm)
```

Apply divided difference operator for `perm` to self.
Linear extension of RCGraph.divdiff_perm.

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.divdiff"></a>

#### divdiff

```python
def divdiff(*seq)
```

Sequential divided difference operators.

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.vertical_coproduct"></a>

#### vertical\_coproduct

```python
def vertical_coproduct()
```

Coproduct of RC graphs, coincides with the coproduct on Schubert polynomials
and induces the mul product.

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.raising_operator"></a>

#### raising\_operator

```python
def raising_operator(index)
```

Linear extension of RCGraph.raising_operator:
Apply raising_operator(index) to every basis RCGraph in self, collect results.
Returns an AltRCGraphRingElement (possibly zero).

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.lowering_operator"></a>

#### lowering\_operator

```python
def lowering_operator(index)
```

Linear extension of RCGraph.lowering_operator.

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.phi"></a>

#### phi

```python
def phi(index)
```

phi(element) := max_{basis rc in supp(element)} phi(rc)
If element is zero, returns 0.

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.epsilon"></a>

#### epsilon

```python
def epsilon(index)
```

epsilon(element) := max_{basis rc in supp(element)} epsilon(rc)

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.crystal_length"></a>

#### crystal\_length

```python
def crystal_length()
```

Use maximum crystal length of basis graphs in support (0 for the zero element).

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.to_highest_weight"></a>

#### to\_highest\_weight

```python
def to_highest_weight()
```

Iteratively raise the element until no further raising is possible.
Returns (highest_weight_element, raise_seq).

Behavior notes:
- This is the natural linear-extension of CrystalGraph.to_highest_weight.
- The returned `highest_weight_element` is an AltRCGraphRingElement.
- raise_seq is the sequence of row indices applied (in order).

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.reverse_raise_seq"></a>

#### reverse\_raise\_seq

```python
def reverse_raise_seq(raise_seq)
```

Apply lowering_operator in reverse order to `raise_seq`.
If the path dies (result is zero), return None (mirrors scalar behavior).

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.crystal_reflection"></a>

#### crystal\_reflection

```python
def crystal_reflection(index)
```

Linear extension of RCGraph.crystal_reflection:
For each basis RCGraph, apply its crystal_reflection(index) and collect results.

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRing"></a>

## AltRCGraphRing Objects

```python
class AltRCGraphRing(SchubertMonomialRing, CrystalGraphRing)
```

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRing.schub"></a>

#### schub

```python
def schub(perm, n=None)
```

Return the AltRCGraphRing element corresponding to the Schubert polynomial
indexed by `perm` in `S_n` (if n is None, n = len(perm) is used).

