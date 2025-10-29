from functools import cache

# import schubmult.rings.free_algebra as fa
from schubmult.schub_lib.perm_lib import Permutation, uncode
from schubmult.symbolic import CoercionFailed, S, sympy_Mul
from schubmult.utils.perm_utils import has_bruhat_descent, mu_A

from .base_schubert_ring import BaseSchubertElement, BaseSchubertRing


def complete_sym_positional_perms_down(orig_perm, p, *k, hack_off=None):
    k = {i - 1 for i in k}
    orig_perm = Permutation(orig_perm)
    total_list = {(orig_perm, 0, 1)}
    up_perm_list = {(orig_perm, 1)}
    # print(f"{orig_perm=}")
    for pp in range(p):
        perm_list = set()
        for up_perm, sign in up_perm_list:
            # pos_list = [i for i in range(k) if up_perm[i] < last]
            rg = [q for q in range(len(up_perm) if hack_off is None else min(len(up_perm),hack_off)) if q not in k and up_perm[q] == orig_perm[q]]
            for j in rg:
                for i in k:
                    a, b = (i, j) if i < j else (j, i)
                    if has_bruhat_descent(up_perm, a, b):
                        # print(f"bruhat ascent on {up_perm=} at {(a,b)=}")
                        new_perm_add = up_perm.swap(a, b)
                        # print(f"{up_perm.inv - orig_perm.inv=}")
                        # print(f"{new_perm_add.inv - orig_perm.inv=}")
                        new_sign = sign if i < j else -sign
                        perm_list.add((new_perm_add, new_sign))
                        total_list.add((new_perm_add, pp + 1, new_sign))
        up_perm_list = perm_list
#     return total_list


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
    # print(f"{bigmu=}, {mu1=}, {mu2=}, {pmu1=}, {pmu2=}, {pmu=} {perm=}, {perm2=}")
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


@cache
def _single_coprod(p, n, T):
    res = T.zero
    for i in range(p + 1):
        res += T.from_dict({((uncode([i]), n), (uncode([p - i]), n)): S.One})
    return res


def _is_code1(perm):
    return perm.inv > 0 and perm.code[0] == perm.inv


class SeparatedDescentsRing(BaseSchubertRing):
    @property
    def args(self):
        return ()

    @property
    def schub_ring(self):
        return self._schub_ring

    # pieri formula for uncode([p]), 1
    def _single_coprod_test(self, p, tensor_elem):
        res = tensor_elem.ring.zero
        for (t1, t2), val in tensor_elem.items():
            for i in range(p + 1):
                # res += T.from_dict({((uncode([i]), n), (uncode([p-i]), n)): S.One})
                telem1 = self.pieri_formula(i, self(*t1))
                telem2 = self.pieri_formula(p - i, self(*t2))
                for perm1, val1 in telem1.items():
                    for perm2, val2 in telem2.items():
                        res += val * val1 * val2 * tensor_elem.ring((perm1, perm2))
        return res

    def pieri_formula(self, p, elem):
        val = self.zero
        for (perm, num_vars), coeff in elem.items():
            lne = len(perm)
            if perm.inv == 0:
                lne = 0
            code_add = p + max(lne, num_vars)
            big_elem = uncode([code_add, *perm.code])
            e_list = complete_sym_positional_perms_down(big_elem, code_add - p, 1)

            for to_add, deg, _ in e_list:
                if deg == code_add - p:
                    if to_add.inv == 0 or max(to_add.descents()) + 1 <= num_vars + 1:
                        val += coeff * self(to_add, num_vars + 1)
        return val

    def __init__(self, ring):
        self._schub_ring = ring
        super().__init__(self._schub_ring.genset, self._schub_ring.coeff_genset)
        self.zero_monom = (self._schub_ring.zero_monom, 0)
        self.dtype = type("SeparatedDescentsRingElement", (SeparatedDescentsRingElement,), {"ring": self})

    def __hash__(self):
        return hash((self._schub_ring, "ARBLE"))

    # @property
    # def rings(self):
    #     return self._rings

    # def domain_new(self, elem1, elem2):
    # @cache
    # def coproduct(self, key):
    #     # R = self @ self
    #     R = fa.FreeAlgebra()
    #     return R.tensor_schub_expand(R.schub_elem(*key).coproduct())

    # @cache
    # def free_element(self, perm, numvars):
    #     return fa.FreeAlgebra().schub_elem(perm, numvars)

    def coproduct_test(self, key):
        T = self @ self
        # if val == self.zero:
        #     return T.zero
        val = self(*key)

        cprd_val = T.zero

        while val != val.ring.zero:
            mx = [k[0].code for k in val.keys() if val[k] != S.Zero]
            mx.sort(reverse=True)
            cd = mx[0]

            mx_key = next(iter([k for k in val.keys() if k[0].code == cd]))
            if len(cd) == 0:
                return cprd_val + T.from_dict({((Permutation([]), mx_key[1]), (Permutation([]), mx_key[1])): val[mx_key] * S.One})
            cd = [*cd]
            fv = cd.pop(0)
            while len(cd) > 1 and cd[-1] == 0:
                cd.pop()
            cf = val[mx_key]
            cprd_val += cf * self._single_coprod_test(fv, self.coproduct_test((uncode(cd), mx_key[1] - 1)))
            val -= cf * self.pieri_formula(fv, self(uncode(cd), mx_key[1] - 1))
        return cprd_val

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
        # return Symbol(f"ASx({k[0]}, {k[1]})", commutative=False)
        return self._schub_ring.printing_term(k, prefix="A")

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
            return [S.Zero]
        return [self[k] if k == self.ring.zero_monom else sympy_Mul(self[k], self.ring.printing_term(k)) for k in self.keys()]
