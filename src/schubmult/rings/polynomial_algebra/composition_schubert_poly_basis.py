from schubmult.combinatorics.permutation import Permutation, uncode
from schubmult.rings.printing import GenericPrintingTerm

from .schubert_poly_basis import SchubertPolyBasis


class CompositionSchubertPolyBasis(SchubertPolyBasis):
    """Wrapper basis for Schubert polynomials indexed by weak compositions.

    Keys are weak compositions interpreted as Lehmer codes. The underlying
    computations are delegated to :class:`SchubertPolyBasis`, while display
    uses the same Schubert printing as the corresponding permutation key.
    """

    def is_key(self, x):
        if isinstance(x, tuple | list) and all(isinstance(a, int) and a >= 0 for a in x):
            return True
        if isinstance(x, tuple | list) and len(x) == 2 and isinstance(x[1], int) and x[1] >= 0:
            return isinstance(x[0], Permutation | tuple | list)
        return False

    def as_key(self, x):
        if isinstance(x, Permutation):
            return tuple(x.trimcode)
        if isinstance(x, tuple | list) and len(x) == 2 and isinstance(x[1], int):
            perm = x[0] if isinstance(x[0], Permutation) else Permutation(x[0])
            return tuple(perm.pad_code(x[1]))
        if not self.is_key(x):
            raise TypeError(f"Invalid composition key: {x!r}")
        return tuple(x)

    def printing_term(self, k):
        return GenericPrintingTerm(f"CompSchub{self.as_key(k)}", "")

    @property
    def zero_monom(self):
        return ()

    @staticmethod
    def _schub_key_to_comp(key):
        perm, length = key
        return tuple(perm.pad_code(length))

    @staticmethod
    def _comp_to_schub_key(comp):
        comp = tuple(comp)
        return (uncode(comp), len(comp))

    def _to_schub_dict(self, dct):
        return {self._comp_to_schub_key(k): v for k, v in dct.items()}

    def _from_schub_dict(self, dct):
        return {self._schub_key_to_comp(k): v for k, v in dct.items()}

    def coproduct(self, key):
        schub_key = self._comp_to_schub_key(key)
        cprd = super().coproduct(schub_key)
        return {((self._schub_key_to_comp(k1), self._schub_key_to_comp(k2))): v for (k1, k2), v in cprd.items()}

    def product(self, key1, key2, coeff=1):
        left = self._comp_to_schub_key(key1)
        right = self._comp_to_schub_key(key2)
        return self._from_schub_dict(super().product(left, right, coeff=coeff))

    def transition(self, other_basis):
        if isinstance(other_basis, CompositionSchubertPolyBasis):
            return lambda x: dict(x)

        return lambda x: self._schubert_transition(other_basis, x)

    def _schubert_transition(self, other_basis, dct):
        schub_dct = self._to_schub_dict(dct)
        transitioned = super().transition(other_basis)(schub_dct)
        return transitioned

    @classmethod
    def dual_basis(cls):
        return SchubertPolyBasis.dual_basis()
