from schubmult.rings.tensor_ring import TensorRing, TensorRingElement

from .crystal_graph import CrystalGraph, CrystalGraphTensor


class CrystalTensorRing(TensorRing):
    def dtype(self):
        elem = CrystalTensorRingElement()
        elem.ring = self
        return elem

def _make_cgt(a, b):
    """
    Create a CrystalGraphTensor for two factors a, b.
    Try common constructor shapes found in repo variants.
    """
    return CrystalGraphTensor(a, b)

def _extract_factors_from_cgt(obj):
    """
    Given a CrystalGraphTensor-like return value, try to extract its factor tuple.
    Prefer .factors attribute, otherwise try to iterate.
    """
    if obj is None:
        return None
    if hasattr(obj, "factors"):
        return tuple(obj.factors)
    try:
        # Some implementations make the object iterable (tuple-like)
        return tuple(obj)
    except Exception:
        return None

class CrystalTensorRingElement(TensorRingElement, CrystalGraph):
    """
    Two-factor tensor ring element (keys are 2-tuples (g1, g2)).
    Delegates two-factor tensor statistics and operators to `CrystalGraphTensor`
    so existing repo tie-breaking is preserved.
    """


    # --- φ / ε : max over support of per-basis two-factor phi/epsilon ---
    def phi(self, index: int) -> int:
        if len(self) == 0:
            return 0
        m = 0
        for (g1, g2) in self.keys():
            try:
                cgt = _make_cgt(g1, g2)
                term = cgt.phi(index)
            except Exception:
                # Defensive fallback
                try:
                    term = max(getattr(g1, "phi", lambda i: 0)(index),
                               getattr(g2, "phi", lambda i: 0)(index))
                except Exception:
                    term = 0
            if term > m:
                m = term
        return m

    def epsilon(self, index: int) -> int:
        if len(self) == 0:
            return 0
        m = 0
        for (g1, g2) in self.keys():
            try:
                cgt = _make_cgt(g1, g2)
                term = cgt.epsilon(index)
            except Exception:
                try:
                    term = max(getattr(g1, "epsilon", lambda i: 0)(index),
                               getattr(g2, "epsilon", lambda i: 0)(index))
                except Exception:
                    term = 0
            if term > m:
                m = term
        return m

    # --- raising/lowering linearized using CrystalGraphTensor ---
    def raising_operator(self, index: int):
        """
        Linear extension: for each basis (g1,g2) create the 2-factor tensor,
        call its raising_operator(index). If non-None, extract new factors
        and add coeff * ring((new_g1, new_g2)) to the result.
        """
        res = self.ring.zero
        for (g1, g2), coeff in self.items():
            try:
                cgt = _make_cgt(g1, g2)
                out = cgt.raising_operator(index)
                new_factors = _extract_factors_from_cgt(out)
                if new_factors:
                    res += coeff * self.ring(tuple(new_factors))
            except Exception:
                # fall back to naive per-factor attempt (should rarely happen)
                try:
                    new_g1 = g1.raising_operator(index)
                    if new_g1 is not None:
                        res += coeff * self.ring((new_g1, g2))
                    else:
                        new_g2 = g2.raising_operator(index)
                        if new_g2 is not None:
                            res += coeff * self.ring((g1, new_g2))
                except Exception:
                    pass
        return res

    def lowering_operator(self, index: int):
        """
        Linear extension of lowering: delegate to CrystalGraphTensor lowering_operator.
        """
        res = self.ring.zero
        for (g1, g2), coeff in self.items():
            try:
                cgt = _make_cgt(g1, g2)
                out = cgt.lowering_operator(index)
                new_factors = _extract_factors_from_cgt(out)
                if new_factors:
                    res += coeff * self.ring(tuple(new_factors))
            except Exception:
                # fallback per-factor attempt
                try:
                    new_g1 = g1.lowering_operator(index)
                    if new_g1 is not None:
                        res += coeff * self.ring((new_g1, g2))
                    else:
                        new_g2 = g2.lowering_operator(index)
                        if new_g2 is not None:
                            res += coeff * self.ring((g1, new_g2))
                except Exception:
                    pass
        return res

    # --- helpers ---
    def crystal_length(self):
        if len(self) == 0:
            return 0
        m = 0
        for (g1, g2) in self.keys():
            l = 0
            try:
                cgt = _make_cgt(g1, g2)
                l = cgt.crystal_length()
            except Exception:
                l = max(getattr(g1, "crystal_length", lambda: len(g1))(),
                        getattr(g2, "crystal_length", lambda: len(g2))())
            if l > m:
                m = l
        return m

    def to_highest_weight(self):
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
