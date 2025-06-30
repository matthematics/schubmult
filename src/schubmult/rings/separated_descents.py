from schubmult.perm_lib import Permutation, uncode
from schubmult.symbolic import S, Symbol, sympy_Mul
from schubmult.utils.perm_utils import add_perm_dict, mu_A, p_trans

from .base_schubert_ring import BaseSchubertElement, BaseSchubertRing


def _sep_desc_mul(perm, perm2, p, q, coeff, ring):
    
    # pmu2 = (~perm2).minimal_dominant_above()
    # mu2 = pmu2.code

    # pmu1 = (~perm).minimal_dominant_above()
    # mu1 = pmu1.code
    # while len(mu2) > 0 and mu2[-1] == 0:
    #     mu2.pop()
    # while len(mu1) > 0 and mu1[-1] == 0:
    #     mu1.pop()

    # if len(mu1) == 0 or mu1[0] < p:
    #     if len(mu1) == 0:
    #         mu1 = list(range(p, 0, -1))
    #     else:
    #         start = mu1[0]
    #         for i in range(p - start):
    #             mu1 = [m + 1 for m in mu1]
    #             mu1 += [1]
    # #mu1 = [max(p,mu1[0])] * q + mu1
    # if len(mu2) == 0:
    #     mu2 = list(range(q, 0, -1))
    
    # if len(mu1) == 0 or mu1[0] < p:
    #     mu1 = [p] + mu1
    # if len(mu2) == 0 or mu2[0] < q:
    #     mu2 = [q] + mu2
    # # if len(mu2) == 0:
    # #     mu2 = [q]
    # # else:
    # #     # start = mu2[0]
    # #     # to_add = q - start
    # #     # mu2 = [a + to_add for a in mu2]
    # #     mu2 = [q] + mu2

    # #mu2 = (max(q - len(mu2), 0)* [start]) + mu2
    # # mu1 = (len(mu2) * [p]) + mu1
    # bigmu = [*mu1]

    # for i in range(len(mu2)):
    #     if i < len(bigmu):
    #         bigmu[i] += mu2[i]
    #     else:
    #         bigmu += [mu2[i]]

    c1 = perm.code
    while len(c1) < p:
        c1 += [0]
    c2 = perm2.code
    while len(c2) < q:
        c2 += [0]
    c = c1 + c2
    n = max([i + c[i] for i in range(len(c))])
    bigmu = list(range(n, 0, -1))


    #print(f"{bigmu=} {p=}, {q=}, {perm=}, {perm2=}, {pmu1=}, {pmu2=}")
    #bigmu = [0] * (p + q)
    # bigmu = []
    # for i in range(p + q, p - 1, -1):
    #     bigmu += [i]
    # bigmu = p_trans(bigmu)

    # if bigmu[0] < p + q:
    #     for i in range(p + q - bigmu[0]):
    #         bigmu = [b+1 for b in bigmu]
            # bigmu += [1]
    mu1 = mu_A(bigmu, list(range(p)))
    mu2 = mu_A(bigmu, list(range(p, len(bigmu))))

    pmu1 = uncode(mu1)
    pmu2 = uncode(mu2)
    pmu = uncode(bigmu)
    #print(f"{bigmu=}, {mu1=}, {mu2=}, {pmu1=}, {pmu2=}, {pmu=} {perm=}, {perm2=}")
    assert (perm * pmu1).inv == pmu1.inv - perm.inv
    assert (perm2 * pmu2).inv == pmu2.inv - perm2.inv
    bingo = ring.from_dict({perm * pmu1: coeff}) * ring.from_dict({perm2 * pmu2: S.One})
    dct = {}
    ipmu = ~pmu
    for w, v in bingo.items():
        if (w * ipmu).inv != ipmu.inv - w.inv:
            continue
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
        ret = self.zero
        for k1, v1 in elem1.items():
            for k2, v2 in elem2.items():
                perm1 = k1[0]
                perm2 = k2[0]

                deg1 = k1[1]
                deg2 = k2[1]

                dct = _sep_desc_mul(perm1, perm2, deg1, deg2, v1 * v2, self.schub_ring)
                dct_update = {}
                for k, v in dct.items():
                    if k.inv > 0 and max(k.descents()) + 1 > deg1 + deg2:
                        continue
                    dct_update[(k, deg1 + deg2)] = v
                ret += self.from_dict(dct_update)
        #return self.from_dict(ret_dict)
        return ret

    def _coerce_add(self, x):  # noqa: ARG002
        return None

    def _coerce_mul(self, x):  # noqa: ARG002
        # have to pull out gens
        # if not isinstance(x.ring, TensorRing):
        #     if set(x.ring.genset) == set(self.genset):
        #         return x.coproduct(*[x.ring.genset.index(v) for v in self.rings[0].genset[1:]])
        return None

    def printing_term(self, k):
        if Permutation.print_as_code:
            return Symbol(f"ASx({k[0].code}, {k[1]})", commutative=False)
        return Symbol(f"ASx({k[0]}, {k[1]})", commutative=False)

    @property
    def one(self):
        return self.from_dict({self.zero_monom: S.One})

    def __call__(self, perm, deg):
        if deg < 0:
            raise ValueError("Degree must be non-negative")
        perm = Permutation(perm)
        if perm.inv > 0:
            deg = max(deg, max(perm.descents()) + 1)
        return self.from_dict({(perm, deg): S.One})


class SeparatedDescentsRingElement(BaseSchubertElement):
    @property
    def free_symbols(self):
        return set()

    def as_ordered_terms(self, *_, **__):
        if len(self.keys()) == 0:
            return [S.Zero]
        return [self[k] if k == self.ring.zero_monom else sympy_Mul(self[k], self.ring.printing_term(k)) for k in self.keys()]
