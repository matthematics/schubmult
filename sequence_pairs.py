from schubmult import *
from schubmult.schub_lib.schub_lib import pull_out_var
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
    mul2 = DSx(perm2*pmu2)
    for perm in elem:
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
        bingo = DSx(perm*pmu1) * mul2
        dct = {}
        ipmu = ~pmu
        for w, v in bingo.items():
            dct[w*ipmu] = v
        res += bingo.ring.from_dict(dct)
    return res

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

pows = [[2], [1], [0, 2]]

res = DSx([])
Permutation.print_as_code = True
for cd in pows:
    res = sep_desc_mul(res, uncode(cd))
    print(res)

        
