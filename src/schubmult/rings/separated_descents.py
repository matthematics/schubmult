from schubmult.schub_lib.permutation import Permutation, uncode
from schubmult.symbolic import CoercionFailed, S, sympify, sympify_sympy, sympy_Mul
from schubmult.utils.perm_utils import mu_A

from .base_schubert_ring import BaseSchubertElement, BaseSchubertRing


def _sep_desc_mul(perm, perm2, p, q, coeff, ring):
    c1 = perm.code
    while len(c1) < p:
        c1 += [0]
    c2 = perm2.code
    while len(c2) < q:
        c2 += [0]
    c = c1 + c2
    if len(c) == 0:
        n = 0
    else:
        n = max([i + c[i] + 1 for i in range(len(c))])
    bigmu = list(range(n, 0, -1))

    mu1 = mu_A(bigmu, list(range(p)))
    mu2 = mu_A(bigmu, list(range(p, len(bigmu))))

    pmu1 = uncode(mu1)
    pmu2 = uncode(mu2)
    pmu = uncode(bigmu)
    assert (perm * pmu1).inv == pmu1.inv - perm.inv
    assert (perm2 * pmu2).inv == pmu2.inv - perm2.inv
    bingo = ring.from_dict({perm * pmu1: coeff}) * ring.from_dict({perm2 * pmu2: S.One})
    dct = {}
    ipmu = ~pmu
    for w, v in bingo.items():
        if (w * ipmu).inv != ipmu.inv - w.inv:
            continue
        # if ring.coeff_genset.label:
        #     v = xreplace_genvars(v, ring.coeff_genset, [-c for c in ring.coeff_genset])
        dct[w * ipmu] = v * (S.NegativeOne ** ((w * ipmu).inv - perm.inv - perm2.inv))
        # xreplace# * (S.NegativeOne ** ((w*ipmu).inv - perm.inv - perm2.inv))
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
        return hash((self._schub_ring, "ARBLE"))

    def mul(self, elem1, elem2):
        # print(f"{elem1=}, {elem2=}")
        try:
            bongus = self.domain_new(elem1)
            return self.from_dict({k: v * bongus for k, v in elem2.items()})
        except CoercionFailed:
            # import traceback
            # traceback.print_exc()
            pass
        try:
            bongus = self.domain_new(elem2)
            return self.from_dict({k: v * bongus for k, v in elem1.items()})
        except CoercionFailed:
            # import traceback
            # traceback.print_exc()
            pass
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
        # return self.from_dict(ret_dict)
        return ret

    def _coerce_add(self, x):
        try:
            x = self.domain_new(x)
            return self.from_dict({self.zero_monom: x})
        except CoercionFailed:
            # import traceback
            # traceback.print_exc()
            return None
        # return None

    def _coerce_mul(self, x):  # noqa: ARG002
        # have to pull out gens
        # if not isinstance(x.ring, TensorRing):
        #     if set(x.ring.genset) == set(self.genset):
        #         return x.coproduct(*[x.ring.genset.index(v) for v in self.rings[0].genset[1:]])
        return None

    def printing_term(self, k):
        # k is a tuple (perm, length)
        from schubmult.rings.abstract_schub_poly import SepDescSchubPoly

        coeff_label = None
        if self.coeff_genset is not NotImplemented and self.coeff_genset is not None:
            coeff_label = self.coeff_genset.label
        return SepDescSchubPoly(k, self.genset.label, coeff_label)

    @property
    def one(self):
        return self.from_dict({self.zero_monom: S.One})

    def _get_min_deg(self, perm):
        if perm.inv > 0:
            return max(perm.descents()) + 1
        return 0

    def new(self, perm, deg=0):
        if not isinstance(perm, Permutation) and not isinstance(perm, list) and not isinstance(perm, tuple):
            elem = self._schub_ring(perm)
            return self.from_dict({(k, self._get_min_deg(k)): v for k, v in elem.items()})
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

    # def coproduct(self):
    #     T = self.ring @ self.ring
    #     res = T.zero
    #     for key, val in self.items():
    #         res += val * self.ring.coproduct(key)
    #     return res

    def coproduct_test(self):
        return self.ring.coproduct_test(self)

    # def free_element(self):
    #     ret = fa.FA([]).ring.zero
    #     for k, v in self.items():
    #         ret += v * self.ring.free_element(*k)
    #     return ret

    def as_ordered_terms(self, *_, **__):
        if len(self.keys()) == 0:
            return [sympify(S.Zero)]
        # Keys are (perm, length) tuples - sort by perm.inv and perm, then by length
        sorted_keys = sorted(self.keys(), key=lambda k: (k[0].inv, tuple(k[0]), k[1]) if k != self.ring.zero_monom else (-1, (), 0))
        return [sympify_sympy(self[k]) if k == self.ring.zero_monom else sympy_Mul(sympify_sympy(self[k]), self.ring.printing_term(k)) for k in sorted_keys]
