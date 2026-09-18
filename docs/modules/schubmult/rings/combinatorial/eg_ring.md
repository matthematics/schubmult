<a id="schubmult.rings.combinatorial.eg_ring"></a>

# schubmult.rings.combinatorial.eg\_ring

`EGRing`: a ring whose basis elements are ``(NilPlactic, length)`` pairs -- the Edelman-Greene
insertion tableau of an RC graph's word together with its row count. ``from_rc_graph`` maps an RC
graph to its EG class.

<a id="schubmult.rings.combinatorial.eg_ring.EGPrintingTerm"></a>

## EGPrintingTerm Objects

```python
class EGPrintingTerm(PrintingTerm)
```

Display symbol for an `EGRing` basis key (prints the key directly).

<a id="schubmult.rings.combinatorial.eg_ring.EGRingElement"></a>

## EGRingElement Objects

```python
class EGRingElement(BaseRingElement)
```

Linear combination of ``(NilPlactic, length)`` basis keys.

<a id="schubmult.rings.combinatorial.eg_ring.EGRing"></a>

## EGRing Objects

```python
class EGRing(BaseRing)
```

The Edelman-Greene tableau ring; see the module docstring.

