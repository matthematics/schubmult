from .abstract_schub_poly import GenericPrintingTerm
from .base_schubert_ring import BaseSchubertElement, BaseSchubertRing
from .free_algebra import FreeAlgebra
from .free_algebra_basis import SchubertBasis, WordBasis
from .rc_graph import RCGraph

FA = FreeAlgebra(WordBasis)
ASx = FreeAlgebra(SchubertBasis)

class DischubertAlgebraElement(BaseSchubertElement):
    pass

class DischubertAlgebra(BaseSchubertRing):

    def __hash__(self):
        return hash((DischubertAlgebra, "pigeon"))

    def printing_term(self, k):
        schubert_part, word_part = k
        return GenericPrintingTerm((schubert_part, word_part), name="ASx⊗FA")

    def __init__(self):
        super().__init__([], [])

    def is_key_valid(self, key):
        return key[0][1] == len(key[1]) and key[0] in FA(*key[1]).change_basis(SchubertBasis)

    def dtype(self):
        elem = DischubertAlgebraElement()
        elem.ring = self
        return elem

    def mul(self, elem1, elem2):
        result = self.zero
        for (schub1, word1), coeff1 in elem1.items():
            for (schub2, word2), coeff2 in elem2.items():
                schub_prod = ASx(*schub1) * ASx(*schub2)
                word_prod = tuple(*word1,*word2)
                for schub_key, schub_coeff in schub_prod.items():
                    new_key = (schub_key, word_prod)
                    result += coeff1 * coeff2 * schub_coeff * self(new_key)
        return result

    def from_dict(self, d):
        elem = DischubertAlgebraElement()
        elem.ring = self
        elem.update({k: v for k, v in d.items() if v != 0})
        return elem

    def coproduct_on_basis(self, k):
        schubert_part, word_part = k
        schub_coprod = ASx(*schubert_part).change_basis(WordBasis).coproduct()
        word_coprod = FA(*word_part).coproduct()
        tring = self @ self
        result = tring.zero
        for (k1, k2), coeff1 in schub_coprod.items():
            for (w1, w2), coeff2 in word_coprod.items():
                elem1 = self.from_dict({(k, w1): v for k, v in FA(*k1).change_basis(SchubertBasis)})
                elem2 = self.from_dict({(k, w2): v for k, v in FA(*k2).change_basis(SchubertBasis)})
                result += coeff1 * coeff2 * tring.ext_multiply(elem1, elem2)
        return result
