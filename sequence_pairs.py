from schubmult import *
from schubmult.schub_lib.schub_lib import pull_out_var
from schubmult.rings import *
from schubmult.symbolic import *
from schubmult.utils.perm_utils import mu_A
from itertools import permutations
import sys

N = int(sys.argv[1])

perms = [Permutation(perm) for perm in list(permutations(list(range(1,N+1))))]

# for perm in perms:
#     if perm.inv == 0:
#         continue
#     num_desc = max(perm.descents()) + 1
#     var_perms = [perm0 for perm0 in perms if len(perm0) <= num_desc]

# max_terms = {}


# for perm in perms:
#     if perm.inv == 0:
#         continue
#     num_terms = len([k for k, v in Sx(perm).as_polynomial().expand().as_coefficients_dict().items() if v != S.Zero])
#     num_vars = max(perm.descents())+1
#     key = (perm.inv, num_vars)
#     max_terms[key] = max(max_terms.get(key, 0),num_terms)

# for (iv, vn), nt in max_terms.items():
#     print(f"{nt}: deg {iv} nv {vn}")


def sep_desc_mul(elem, perm2):
    pmu2 = (~perm2).minimal_dominant_above()
    mu2 = pmu2.code
    res = S.Zero
    mul2 = elem.ring.from_dict({perm2*pmu2: S.One})
    for perm, v2 in elem.items():
        pmu1 = (~perm).minimal_dominant_above()
        mu1 = pmu1.code
        if len(mu1) > 0:
            mu1 = (len(mu2)*[mu1[0]])+mu1
        bigmu = [*mu1]

        for i in range(len(mu2)):
            if i<len(bigmu):
                bigmu[i] += mu2[i]
            else:
                bigmu += [mu2[i]]
        pmu1 = uncode(mu1)
        pmu = uncode(bigmu)
        bingo = mul2.ring.from_dict({perm*pmu1: v2}) * mul2
        dct = {}
        ipmu = ~pmu
        for w, v in bingo.items():
            dct[w*ipmu] = v
        res += bingo.ring.from_dict(dct)
    return res

def h_coprod(p, k):
    ret = []
    for i in range(p + 1):
        if i == 0:
            ret += [[[], ([0] * (k-1)) + [ p -i ]]]
        elif i == p:
            ret += [[([0] * (k-1)) + [ p ], []]]
        else:
            ret += [[([0] * (k-1)) + [ i], ([0] * (k-1)) + [p - i]]]
    return ret
    # for perm2 in perms:
    #     mu1 = (~perm).minimal_dominant_above().code
    #     mu2 = (~perm2).minimal_dominant_above().code
    #     if len(mu1) > 0:
    #         mu1 = (len(mu2)*[mu1[0]])+mu1
    #     bigmu = [*mu1]

    #     for i in range(len(mu2)):
    #         if i<len(bigmu):
    #             bigmu[i] += mu2[i]
    #         else:
    #             bigmu += [mu2[i]]
    #     mu1 = uncode(mu1)
    #     mu2 = uncode(mu2)
    #     mu = uncode(bigmu)
    #     bingo = DSx(perm*mu1) * DSx(perm2*mu2)
    #     dct = {}
    #     for w, v in bingo.items():
    #         dct[w*(~mu)] = v
    #     print(f"{perm}*{perm2} = {bingo.ring.from_dict(dct)}")

  # pows = [[0, 0, 1], [3], [5]]
pows = [[1, 3], [2, 1]]


res = Sx([])

RR = TensorRing(res.ring, res.ring)

elems = []


for pork in pows:
    dct = {}
    for panko, fong in h_coprod(*pork):
        dct[(uncode(panko), uncode(fong))] = S.One
    elems += [RR.from_dict(dct)]

Permutation.print_as_code = True
res = RR.one
from schubmult.rings._mul_utils import _tensor_product_of_dicts, add_perm_dict

for elem in elems:
    ret = S.Zero
    for k1, v1 in res.items():
        for k2, v2 in elem.items():
            dct = sep_desc_mul(RR.rings[0].from_dict({k1[0]: v1 * v2}), k2[0])
            dct2 = sep_desc_mul(RR.rings[1].from_dict({k1[1]: 1}),k2[1])
            dct = _tensor_product_of_dicts(dct, dct2)
            ret += RR.from_dict(dct)
    res = ret
    print(res)
# for cd in pows:
#     res = sep_desc_mul(res, uncode(cd))
#     print(res.expand(deep=False))

        
