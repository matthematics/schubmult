# src/schubmult/rings/crystal_graph_ring.py



## class CrystalGraphRing(BaseSchubertRing)

Ring whose basis elements are CrystalGraph-like objects.

We deliberately do not special-case tensor objects here: CrystalGraphTensor
implements the same CrystalGraph API and will be handled by polymorphism.

## _ensure_cg(obj)

Ensure obj behaves like a CrystalGraph. If it's already a CrystalGraph,
return it. Otherwise return obj unchanged and rely on duck-typing by callers.

## class CrystalGraphRingElement(BaseSchubertElement, CrystalGraph)

Element of the CrystalGraphRing.

Keys are arbitrary objects that implement the CrystalGraph API (including
CrystalGraphTensor). All crystal operators / statistics are lifted linearly
by delegating to the underlying key's methods.

