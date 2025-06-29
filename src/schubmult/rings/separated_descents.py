from schubmult.perm_lib import Permutation, uncode
from schubmult.symbolic import S, Symbol, sympy_Mul
from schubmult.utils.perm_utils import add_perm_dict, mu_A

from .base_schubert_ring import BaseSchubertElement, BaseSchubertRing


def _sep_desc_mul(perm, perm2, p, coeff, ring):
    pmu2 = (~perm2).minimal_dominant_above()
    mu2 = pmu2.code
    mul2 = ring.from_dict({perm2 * pmu2: S.One})
    pmu1 = (~perm).minimal_dominant_above()
    mu1 = pmu1.code
    if len(mu1) > 0:
        mu1 = (len(mu2) * [p]) + mu1
    bigmu = [*mu1]

    for i in range(len(mu2)):
        if i < len(bigmu):
            bigmu[i] += mu2[i]
        else:
            bigmu += [mu2[i]]
    mu1 = mu_A(bigmu, list(range(p)))
    mu2 = mu_A(bigmu, list(range(p, len(bigmu))))
    pmu1 = uncode(mu1)
    pmu = uncode(bigmu)
    bingo = ring.from_dict({perm * pmu1: coeff}) * mul2
    dct = {}
    ipmu = ~pmu
    for w, v in bingo.items():
        dct[w * ipmu] = v
    return dct


class SeparatedDescentsRing(BaseSchubertRing):
    @property
    def args(self):
        return ()

    @property
    def schub_ring(self):
        return self._schub_ring

    def __init__(self, ring):
        self._schub_ring = ring
        super().__init__(self._schub_ring.genset, self._schub_ring.coeff_genset)
        self.zero_monom = (self._schub_ring.zero_monom, 0)
        self.dtype = type("SeparatedDescentsRingElement", (SeparatedDescentsRingElement,), {"ring": self})

    def __hash__(self):
        return hash(self._rings)

    # @property
    # def rings(self):
    #     return self._rings

    def mul(self, elem1, elem2):
        ret_dict = {}
        for k1, v1 in elem1.items():
            for k2, v2 in elem2.items():
                perm1 = k1[0]
                perm2 = k2[0]

                deg1 = k1[1]
                deg2 = k2[1]

                dct = _sep_desc_mul(perm1, perm2, deg1, v1 * v2, self.schub_ring)
                to_add = {(k, deg1 + deg2): v for k, v in dct.items()}
                ret_dict = add_perm_dict(ret_dict, to_add)
        return self.from_dict(ret_dict)

    def _coerce_add(self, x):  # noqa: ARG002
        return None

    def _coerce_mul(self, x):  # noqa: ARG002
        # have to pull out gens
        # if not isinstance(x.ring, TensorRing):
        #     if set(x.ring.genset) == set(self.genset):
        #         return x.coproduct(*[x.ring.genset.index(v) for v in self.rings[0].genset[1:]])
        return None

    def printing_term(self, k):
        return Symbol(f"ASx({k[0]}, {k[1]})", commutative=False)

    @property
    def one(self):
        return self.from_dict({self.zero_monom: S.One})

    def __call__(self, perm, deg):
        return self.from_dict({(Permutation(perm), deg): S.One})


class SeparatedDescentsRingElement(BaseSchubertElement):
    @property
    def free_symbols(self):
        return set()

    def as_ordered_terms(self, *_, **__):
        if len(self.keys()) == 0:
            return [S.Zero]
        return [self[k] if k == self.ring.zero_monom else sympy_Mul(self[k], self.ring.printing_term(k)) for k in self.keys()]
