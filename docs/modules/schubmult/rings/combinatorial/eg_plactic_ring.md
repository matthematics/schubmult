<a id="schubmult.rings.combinatorial.eg_plactic_ring"></a>

# schubmult.rings.combinatorial.eg\_plactic\_ring

`EGPlacticRing`: a `CrystalGraphRing` on pairs ``((NilPlactic, length), Plactic)`` -- an RC
graph's Edelman-Greene insertion tableau together with its plactic recording tableau -- with
conversion to and from `RCGraphRing`.

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticPrintingTerm"></a>

## EGPlacticPrintingTerm Objects

```python
class EGPlacticPrintingTerm(PrintingTerm)
```

Display symbol for an `EGPlacticRing` basis key.

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticRingElement"></a>

## EGPlacticRingElement Objects

```python
class EGPlacticRingElement(CrystalGraphRingElement)
```

EGPlacticRing elements are linear combinations of ``((NilPlactic, length), Plactic)`` basis keys.

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticRingElement.__mod__"></a>

#### \_\_mod\_\_

```python
def __mod__(other)
```

Polynomial product: self % other.
Currently only defined when `other` is a dominant RC graph.

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticRingElement.raising_operator"></a>

#### raising\_operator

```python
def raising_operator(index)
```

Linear extension of RCGraph.raising_operator:
Apply raising_operator(index) to every basis RCGraph in self, collect results.
Returns an EGPlacticRingElement (possibly zero).

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticRingElement.lowering_operator"></a>

#### lowering\_operator

```python
def lowering_operator(index)
```

Linear extension of RCGraph.lowering_operator.

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticRingElement.phi"></a>

#### phi

```python
def phi(index)
```

phi(element) := max_{basis rc in supp(element)} phi(rc)
If element is zero, returns 0.

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticRingElement.epsilon"></a>

#### epsilon

```python
def epsilon(index)
```

epsilon(element) := max_{basis rc in supp(element)} epsilon(rc)

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticRingElement.crystal_length"></a>

#### crystal\_length

```python
def crystal_length()
```

Use maximum crystal length of basis graphs in support (0 for the zero element).

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticRingElement.to_highest_weight"></a>

#### to\_highest\_weight

```python
def to_highest_weight()
```

Iteratively raise the element until no further raising is possible.
Returns (highest_weight_element, raise_seq).

Behavior notes:
- This is the natural linear-extension of CrystalGraph.to_highest_weight.
- The returned `highest_weight_element` is an EGPlacticRingElement.
- raise_seq is the sequence of row indices applied (in order).

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticRingElement.to_lowest_weight"></a>

#### to\_lowest\_weight

```python
def to_lowest_weight()
```

Iteratively raise the element until no further raising is possible.
Returns (highest_weight_element, raise_seq).

Behavior notes:
- This is the natural linear-extension of CrystalGraph.to_highest_weight.
- The returned `highest_weight_element` is an EGPlacticRingElement.
- raise_seq is the sequence of row indices applied (in order).

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticRingElement.reverse_raise_seq"></a>

#### reverse\_raise\_seq

```python
def reverse_raise_seq(raise_seq)
```

Apply lowering_operator in reverse order to `raise_seq`.
If the path dies (result is zero), return None (mirrors scalar behavior).

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticRing"></a>

## EGPlacticRing Objects

```python
class EGPlacticRing(CrystalGraphRing)
```

The EG-plactic ring; see the module docstring.

