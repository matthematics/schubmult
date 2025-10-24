from typing import Any, Optional, Tuple

from .base_schubert_ring import BaseSchubertElement, BaseSchubertRing
from .crystal_graph import CrystalGraph

# Prefer the explicit BaseSchubertRing if available; fall back to TensorRing to remain safe.

class CrystalGraphRing(BaseSchubertRing):
    """
    Ring whose basis elements are CrystalGraph-like objects.

    We deliberately do not special-case tensor objects here: CrystalGraphTensor
    implements the same CrystalGraph API and will be handled by polymorphism.
    """

    def dtype(self):
        elem = CrystalGraphRingElement()
        elem.ring = self
        return elem


def _ensure_cg(obj: Any) -> Any:
    """
    Ensure obj behaves like a CrystalGraph. If it's already a CrystalGraph,
    return it. Otherwise return obj unchanged and rely on duck-typing by callers.
    """
    if obj is None:
        return None
    if isinstance(obj, CrystalGraph):
        return obj
    return obj


class CrystalGraphRingElement(BaseSchubertElement, CrystalGraph):
    """
    Element of the CrystalGraphRing.

    Keys are arbitrary objects that implement the CrystalGraph API (including
    CrystalGraphTensor). All crystal operators / statistics are lifted linearly
    by delegating to the underlying key's methods.
    """

    def phi(self, index: int) -> int:
        if len(self) == 0:
            return 0
        m = 0
        for key in self.keys():
            try:
                cg = _ensure_cg(key)
                term = cg.phi(index)
            except Exception:
                try:
                    term = getattr(key, "phi", lambda i: 0)(index)
                except Exception:
                    term = 0
            if term > m:
                m = term
        return m

    def epsilon(self, index: int) -> int:
        if len(self) == 0:
            return 0
        m = 0
        for key in self.keys():
            try:
                cg = _ensure_cg(key)
                term = cg.epsilon(index)
            except Exception:
                try:
                    term = getattr(key, "epsilon", lambda i: 0)(index)
                except Exception:
                    term = 0
            if term > m:
                m = term
        return m

    def raising_operator(self, index: int):
        """
        Linearized raising operator: delegate to each key's raising_operator
        and collect results in the ring.
        """
        res = self.ring.zero
        for key, coeff in self.items():
            try:
                cg = _ensure_cg(key)
                out = cg.raising_operator(index)
                if out is not None:
                    res += coeff * self.ring(out)
            except Exception:
                # Defensive fallback: try attribute-based per-key operator
                try:
                    new_key = getattr(key, "raising_operator", lambda i: None)(index)
                    if new_key is not None:
                        res += coeff * self.ring(new_key)
                except Exception:
                    pass
        return res

    def lowering_operator(self, index: int):
        """
        Linearized lowering operator: delegate to each key's lowering_operator.
        """
        res = self.ring.zero
        for key, coeff in self.items():
            try:
                cg = _ensure_cg(key)
                out = cg.lowering_operator(index)
                if out is not None:
                    res += coeff * self.ring(out)
            except Exception:
                try:
                    new_key = getattr(key, "lowering_operator", lambda i: None)(index)
                    if new_key is not None:
                        res += coeff * self.ring(new_key)
                except Exception:
                    pass
        return res

    def crystal_length(self) -> int:
        if len(self) == 0:
            return 0
        m = 0
        for key in self.keys():
            try:
                cg = _ensure_cg(key)
                l = cg.crystal_length()
            except Exception:
                try:
                    l = getattr(key, "crystal_length", lambda: 0)()
                except Exception:
                    l = 0
            if l > m:
                m = l
        return m

    def to_highest_weight(self) -> Tuple["CrystalGraphRingElement", Tuple[int, ...]]:
        elem = self
        seq = []
        found = True
        while found:
            found = False
            for row in range(1, elem.crystal_length()):
                new_e = elem.raising_operator(row)
                if new_e is not None and len(new_e) > 0 and new_e != elem:
                    elem = new_e
                    seq.append(row)
                    found = True
                    break
        return elem, tuple(seq)
