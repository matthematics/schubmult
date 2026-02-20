from functools import cache

import schubmult.mult.double as yz
import schubmult.mult.single as py
import schubmult.rings.schubert.quantum_schubert_ring as qsr
import schubmult.utils.schub_lib as schub_lib
from schubmult.combinatorial_reps.permutation import Permutation
from schubmult.symbolic import S
from schubmult.symbolic.poly.poly_lib import xreplace_genvars
from schubmult.symbolic.poly.variables import GeneratingSet, GeneratingSet_base, poly_genset
from schubmult.symmetric_polynomials import ElemSym

from .double_schubert_ring import DoubleSchubertElement, DoubleSchubertRing, DSx, ElemDoubleSchubertRing

__all__ = [
    "DSx",
    "DoubleSchubertElement",
    "DoubleSchubertRing",
    "ElemDoubleSchubertRing",
    "SingleSchubertRing",
    "Sx",
]


class SingleSchubertRing(DoubleSchubertRing):
    def __init__(self, genset):
        super().__init__(genset, poly_genset(0))
        self.dtype = type("DoubleSchubertElement", (DoubleSchubertElement,), {"ring": self})

    def __str__(self):
        return f"Schubert polynomial ring in {self.genset.label}"

    def __hash__(self):
        return hash((self.genset, self.coeff_genset, "SBS"))

    def _coerce_mul(self, other):
        """Coerce a basis schubert algebra element so it can be multiplied

        Args:
            other (_type_): _description_

        Returns:
            _type_: _description_
        """
        if type(other.ring) is type(self):
            if self.genset == other.ring.genset:
                return other
        if type(other.ring) is DoubleSchubertRing:
            if self.genset == other.ring.genset:
                return other
        if isinstance(other.ring, qsr.QuantumDoubleSchubertRing):
            return self._coerce_mul(other.as_classical())
        return None

    @cache
    def cached_product(self, u, v, basis2):
        if self == basis2:
            return py.schubmult_py({u: S.One}, v)
        return {k: xreplace_genvars(x, poly_genset(0), basis2.coeff_genset if basis2.coeff_genset else poly_genset(0)) for k, x in yz.schubmult_double_pair_generic(u, v).items()}

    @cache
    def cached_positive_product(self, u, v, basis2):
        return self.cached_product(u, v, basis2)

    def single_variable(self, elem, varnum):
        ret = self.zero
        for u, v in elem.items():
            new_perms = schub_lib.elem_sym_positional_perms(u, 1, varnum)
            for perm, udiff, sign in new_perms:
                if udiff == 1:
                    ret += self.domain_new(sign * v) * self.new(perm)
        return ret

    def new(self, x):
        genset = self.genset
        if not genset:
            genset = self.genset
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        if isinstance(x, list) or isinstance(x, tuple):
            elem = self.from_dict({Permutation(x): self.domain.one})
        elif isinstance(x, Permutation):
            elem = self.from_dict({x: S.One})
        else:
            elem = self.from_expr(x)
        return elem

    @property
    def elem_func(self):
        return ElemSym


Sx = SingleSchubertRing(GeneratingSet("x"))
