import sys

from schubmult import Permutation, schub_coprod_py, schubmult_py, uncode
from schubmult import trimcode
from schubmult import S
from schubmult.utils.argparse import schub_argparse
from schubmult.utils.perm_utils import mu_A, p_trans


def is_grass(perm):
    return len(perm.descents()) <= 1

def is_k_separated(perm1, perm2):
    if perm1.inv == 0 or perm2.inv == 0:
        return None
    mn = min(perm2.descents())
    mx = max(perm2.descents())
    for d in perm1.descents():
        if (d > mn and d <= mx) or (d>=mn and d<mx):
            return None
    return mn, mx

def sum_lists(L1, L2):
    L3 = [0] * (max(len(L1), len(L2)))
    if len(L1) > len(L2):
        L1, L2 = L2, L1
    for i in range(len(L1)):
        L3[i] = L1[i] + L2[i]
    for i in range(len(L1), len(L2)):
        L3[i] = L2[i]
    return L3

def generate_AB(mu_u, mu_v):
    t_mu_u_t = (~mu_u).code
    t_mu_v_t = (~mu_v).code
    t_mu_u = (mu_u).code
    t_mu_v = (mu_v).code
    t_mu_uv = sum_lists(t_mu_u, t_mu_v)
    t_mu_uv_t = p_trans(t_mu_uv)

    # t_mu_uv = [t_mu_u[i] + t_mu_v[i] for i in range(len(t_mu_u))]
    # t_mu_uv_t = p_trans(t_mu_uv)
    A = []
    B = []
    indexA = 0

    while len(t_mu_u_t) > 0 and t_mu_u_t[-1] == 0:
        t_mu_u_t.pop()

    while len(t_mu_v_t) > 0 and t_mu_v_t[-1] == 0:
        t_mu_v_t.pop()

    while len(t_mu_uv_t) > 0 and t_mu_uv_t[-1] == 0:
        t_mu_uv_t.pop()

    for index in range(len(t_mu_uv_t)):
        if indexA < len(t_mu_u_t) and t_mu_uv_t[index] == t_mu_u_t[indexA]:
            A += [index]
            indexA += 1
        else:
            B += [index]

    # mu_w_A = uncode(mu_A(mu_w., A))
    # mu_w_B = uncode(mu_A(mu_w.code, B))
    return A, B

if __name__ == "__main__":
    N = int(sys.argv[1])
    perms = Permutation.all_permutations(N)

    for perm1 in perms:
        for perm2 in perms:
            interv = is_k_separated(perm1, perm2)
            if not interv:
                continue
            mn8, mx8 = interv
            mu_u = ~((~perm1).minimal_dominant_above())
            mu_v = ~((~perm2).minimal_dominant_above())

            mu_w = uncode(sum_lists((mu_u).code, (mu_v).code))
            coeff_dict = schubmult_py({perm1: S.One}, perm2)

            #A = list(range(mn+2,mx+2))
            A, B = generate_AB(mu_v, mu_u)
            # if all(a < b for a in A for b in B):
            #     continue
            # if all(a > b for a in A for b in B):
            #     continue
            # end = max(max(A),max(B))+1
            # B = [*B,*list(range(end, end+5))]
            #A_left = list(range(mn,mx+1))
            for w, coeff in coeff_dict.items():
                w_prime = w * (~mu_w)
                u_prime = perm1 * (~mu_u)
                v_prime = perm2 * (~mu_v) # this is A
                mn, mx = min(A)-1, max(A)
                if len(w_prime.descents()) > 3:
                    w_prime_reduced, ddiff_perm = w_prime.parabolic_reduce(*list({mn, mx}))
                    from bisect import bisect_left
                    while ddiff_perm.inv > 0:
                        desc = next(iter(ddiff_perm.descents()))
                        ddiff_perm = ddiff_perm.swap(desc, desc+1)
                        if desc in A:
                            desc2 = bisect_left(A, desc)
                            v_prime = v_prime.swap(desc2, desc2 + 1)
                        elif desc in B:
                            desc2 = bisect_left(B, desc)
                            u_prime = u_prime.swap(desc2, desc2 + 1)
                        else:
                            raise ValueError(f"{desc=} not in {A=} or {B=}")
                        # if desc < mn:
                        #     u_prime = u_prime.swap(desc, desc+1)
                        # elif desc > mx:
                        #     desc2 = desc - mx
                        #     u_prime = u_prime.swap(desc2, desc2+1)
                        # else:
                        #     desc2 = desc - mn
                        #     v_prime = v_prime.swap(desc2, desc2 + 1)
                    assert len(w_prime_reduced.descents()) <= 3
                    new_mu_w = (~w_prime_reduced).minimal_dominant_above()
                    top_coeff = w_prime_reduced * (new_mu_w)
                    assert len(top_coeff.descents()) <= 3
                    #B = [*list(range(1,mn + 2)), *list(range(mx+2, max(w_prime_reduced.descents())+1))]
                    new_mu_u = uncode(mu_A(new_mu_w.code,B))
                    new_mu_v = uncode(mu_A(new_mu_w.code,A))
                    u_coeff = u_prime * new_mu_u
                    v_coeff = v_prime * new_mu_v

                    assert len(u_coeff.descents()) <= 3
                    assert len(v_coeff.descents()) <= 3
                    try:
                        assert schubmult_py({u_coeff: S.One}, v_coeff).get(top_coeff, S.Zero) == coeff
                    except AssertionError:
                        raise


    # argv = sys.argv
    # args, formatter = schub_argparse(
    #         "schubmult_py",
    #         "Compute products of ordinary Schubert polynomials with Schur polys",
    #         argv=argv[1:],
    #     )
    # perms = args.perms
    # if len(perms) > 2:
    #     # print("Dunna do dees too long arghs")
    #     sys.exit(1)

    # to_mul = uncode(perms[0])
    # grass_perm = uncode(perms[1])
    # assert is_grass(grass_perm)

    # cd = trimcode(grass_perm)

    # # c_{u\mu_A, v\mu_B}^{w\mu}
    # # grass_perm = vmub
    # # ~v = mub(~grass_perm)
    # # v = grass_perm~mub
    # muB = ~()

    # # indices = [*list(range(1,1+cd.count(0))), *list(range(len(cd)+1,len(cd)+to_mul.inv))]


    # # kcd = [indices[i] - i - 1 for i in range(len(indices))] + [n + 1 - k for i in range(k, n)]
    # # max_required = max([kcd[i] + i for i in range(len(kcd))])
    # # kcd2 = kcd + [0 for i in range(len(kcd), max_required)] + [0]

    # # kperm = ~(uncode(kcd2))




    # # def schub_coprod_py(perm, indices):
    # # mperm = perm
    # # indices = sorted(indices)
    # # ret_dict = {}
    # # k = len(indices)
    # # n = len(mperm)
    # # kcd = [indices[i] - i - 1 for i in range(len(indices))] + [n + 1 - k for i in range(k, n)]
    # # max_required = max([kcd[i] + i for i in range(len(kcd))])
    # # kcd2 = kcd + [0 for i in range(len(kcd), max_required)] + [0]
    # # N = len(kcd)
    # # kperm = ~(uncode(kcd2))
    # # coeff_dict = {kperm: 1}
    # # logger.debug(f"{kperm.code}*{mperm.code}")
    # # coeff_dict = schubmult_py(coeff_dict, mperm)

    # # inv_kperm = inv(kperm)
    # # inverse_kperm = ~kperm
    # # # total_sum = 0
    # # for perm, val in coeff_dict.items():
    # #     if val == 0:
    # #         continue
    # #     # pperm = [*perm]
    # #     downperm = perm * inverse_kperm
    # #     if inv(downperm) == inv(perm) - inv_kperm:
    # #         flag = True
    # #         for i in range(N):
    # #             if downperm[i] > N:
    # #                 flag = False
    # #                 break
    # #         if not flag:
    # #             continue
    # #         firstperm = permtrim(list(downperm[0:N]))
    # #         secondperm = permtrim([downperm[i] - N for i in range(N, len(downperm))])
    # #         ret_dict[(firstperm, secondperm)] = val
    # # return ret_dict
