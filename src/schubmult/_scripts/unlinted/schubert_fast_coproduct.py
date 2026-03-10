from schubmult import *

from schubmult.utils.perm_utils import has_bruhat_descent
from schubmult.utils.schub_lib import compute_vpathdicts, elem_sym_perms

# def bruhat_interval(perm, max_code_length=None):
#     return [p for p in Permutation.all_permutations(len(perm)) if p.bruhat_leq(perm) and (max_code_length is None or len(p.trimcode) <= max_code_length)]

def fa_coproduct(perm, length):
    
    m = len(perm)
    w0 = Permutation.w0(m)
    v = w0 * perm
    result = ASx.zero @ ASx.zero
    
    vn1 = ~v
    th = vn1.theta()
    # if len(th) == 0 or th[0] == 0:
    #     return ASx(perm, length) @ ASx(perm, length)
    mu = uncode(th)
    vmu = v * mu
    inv_vmu = vmu.inv
    inv_mu = mu.inv
    while len(th) > 0 and th[-1] == 0:
        th.pop()
    if len(th) == 0:
        vpathdicts = []
    else:
        vpathdicts = compute_vpathdicts(th, vmu)
    mx_th = [0 for i in range(len(th))]
    for index in range(len(th)):
        for vp in vpathdicts[index]:
            for v2, vdiff, s in vpathdicts[index][vp]:
                mx_th[index] = max(mx_th[index], th[index] - vdiff)

    visited = {perm}
    perm_list = [perm]
    while perm_list:
        u = perm_list.pop()
        for i in range(len(u)-1):
            for j in range(i+1, len(u)):
                if has_bruhat_descent(u, i, j):
                    new_u = u.swap(i, j)
                    if new_u not in visited and len(new_u.trimcode) <= length:
                        visited.add(new_u)
                        perm_list.append(new_u)
    #visited = dict.fromkeys(visited, 1)
        #u = perm_list.pop()
    result = sum([ASx(u, length) @ ASx([], length) for u in visited])
    for u in visited:
            # the_prod = Sx(p) * Sx(prod_perm)
            # for q, coeff in the_prod.items():
            #     if len(q) > m:
            #         continue
            #     the_term = w0 * q
            #     if len(the_term.trimcode) > trmcode:
            #         continue
            #     result += coeff * ASx(p, trmcode) @ ASx(the_term, trmcode)
            # for i in range(len(p)-1):
            #     for j in range(i+1, len(p)):
            #         if has_bruhat_descent(p, i, j):
            #             new_p = p.swap(i, j)
            #             if new_p not in visited:
            #                 visited.add(new_p)
            #                 perm_list.append(new_p)
        if v.inv != 0:
            inv_u = u.inv
            vpathsums = {Permutation(u): {Permutation([]): 1}}

            for index in range(len(th)):
                newpathsums = {}
                for up in vpathsums:
                    inv_up = up.inv
                    newperms = elem_sym_perms(
                        up,
                        min(mx_th[index], inv_mu - inv_vmu - (inv_up - inv_u)),
                        th[index],
                    )
                    # print(f"{up=}")
                    for vp in vpathsums[up]:
                        # print(f"{vp=} {type(vp)=} {hash(vp)=}")
                        # print(f"{vpathsums[up]=} {vpathdicts[index]=}")
                        sumval = vpathsums[up][vp]
                        if sumval == 0:
                            continue
                        for v2, vdiff, s in vpathdicts[index][vp]:
                            addsumval = s * sumval
                            for up2, udiff in newperms:
                                if vdiff + udiff == th[index]:
                                    if up2 not in newpathsums:
                                        newpathsums[up2] = {}
                                    newpathsums[up2][v2] = newpathsums[up2].get(v2, 0) + addsumval
                vpathsums = newpathsums
            toget = vmu
            coeff_dict = {ep: vpathsums[ep].get(toget, 0) for ep in vpathsums if len(ep) <= m}
        else:
            coeff_dict = {u: 1}
        for ep in coeff_dict:
            coeff = coeff_dict[ep]
            if coeff != 0:
                the_term = w0 * ep
                result += coeff * ASx(u, length) @ ASx(the_term, length)
        
    return result


if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)
    
    for perm in perms:
        result = ASx.zero @ ASx.zero
        trmcode = len(perm.trimcode) + 2
        result = fa_coproduct(perm, trmcode)
        #m = len(perm)
        #w0 = Permutation.w0(m)
        #prod_perm = w0 * perm
        # interval = bruhat_interval(perm, trmcode)
        # interval_by_length = {}
        # for p in interval:
        #     interval_by_length[p.inv] = interval_by_length.get(p.inv, []) + [p]
        # visited = set()
        # perm_list = [perm]
        # while perm_list:
        #     p = perm_list.pop()
        #     the_prod = Sx(p) * Sx(prod_perm)
        #     for q, coeff in the_prod.items():
        #         if len(q) > m:
        #             continue
        #         the_term = w0 * q
        #         if len(the_term.trimcode) > trmcode:
        #             continue
        #         result += coeff * ASx(p, trmcode) @ ASx(the_term, trmcode)
        #     for i in range(len(p)-1):
        #         for j in range(i+1, len(p)):
        #             if has_bruhat_descent(p, i, j):
        #                 new_p = p.swap(i, j)
        #                 if new_p not in visited:
        #                     visited.add(new_p)
        #                     perm_list.append(new_p)
        assert result.almosteq(ASx(perm, trmcode).coproduct()), f"Failed for {perm} with coproduct \n{result}\nvermin\n{ASx(perm, trmcode).coproduct()}"
        print(f"Passed for {perm}")