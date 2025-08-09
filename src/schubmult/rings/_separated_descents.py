from functools import cache

# import schubmult.rings.free_algebra as fa
from schubmult.perm_lib import uncode
from schubmult.symbolic import S
from schubmult.utils.perm_utils import mu_A


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
