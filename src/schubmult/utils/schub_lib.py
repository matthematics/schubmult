from functools import cache

import schubmult.rings.variables as spl
from schubmult.schub_lib.perm_lib import (
    Permutation,
    code,
    inv,
    inverse,
    one_dominates,
    permtrim,
    theta,
    uncode,
)
from schubmult.symbolic import S, prod
from schubmult.utils.perm_utils import (
    count_bruhat,
    get_cycles,
    getpermval,
    has_bruhat_ascent,
    has_bruhat_descent,
    mu_A,
    omega,
    p_trans,
    sg,
)

q_var = spl.GeneratingSet("q")


# def double_elem_sym_q(u, p1, p2, k, q_var=q_var):
#     ret_list = {}
#     perms1 = elem_sym_perms_q(u, p1, k, q_var)
#     iu = inverse(u)
#     for perm1, udiff1, mul_val1 in perms1:
#         perms2 = elem_sym_perms_q(perm1, p2, k, q_var)
#         cycles1 = get_cycles(tuple(permtrim(mulperm(iu, [*perm1]))))
#         cycles1_dict = {}
#         for c in cycles1:
#             if c[-1] not in cycles1_dict:
#                 cycles1_dict[c[-1]] = []
#             cycles1_dict[c[-1]] += [set(c)]
#         ip1 = inverse(perm1)
#         for perm2, udiff2, mul_val2 in perms2:
#             cycles2 = get_cycles(tuple(permtrim(mulperm(ip1, [*perm2]))))
#             good = True
#             for i in range(len(cycles2)):
#                 c2 = cycles2[i]
#                 if c2[-1] not in cycles1_dict:
#                     continue
#                 for c1_s in cycles1_dict[c2[-1]]:
#                     for a in range(len(c2) - 2, -1, -1):
#                         if c2[a] in c1_s:
#                             good = False
#                             break
#                     if not good:
#                         break
#                 if not good:
#                     break

#             if good:
#                 # print(f"{(perm1, udiff1, mul_val1)=}")
#                 if (perm1, udiff1, mul_val1) not in ret_list:
#                     ret_list[(perm1, udiff1, mul_val1)] = []
#                 ret_list[(perm1, udiff1, mul_val1)] += [(perm2, udiff2, mul_val2)]
#     return ret_list


def double_elem_sym_q(u, p1, p2, k, q_var=q_var):
    ret_list = {}
    perms1 = elem_sym_perms_q(u, p1, k, q_var)
    iu = ~u
    for perm1, udiff1, mul_val1 in perms1:
        perms2 = elem_sym_perms_q(perm1, p2, k, q_var)
        cycles1 = get_cycles(iu * perm1)
        cycles1_dict = {}
        for c in cycles1:
            if c[-1] not in cycles1_dict:
                cycles1_dict[c[-1]] = []
            cycles1_dict[c[-1]] += [set(c)]
        ip1 = inverse(perm1)
        for perm2, udiff2, mul_val2 in perms2:
            cycles2 = get_cycles(ip1 * perm2)
            good = True
            for i in range(len(cycles2)):
                c2 = cycles2[i]
                if c2[-1] not in cycles1_dict:
                    continue
                for c1_s in cycles1_dict[c2[-1]]:
                    for a in range(len(c2) - 2, -1, -1):
                        if c2[a] in c1_s:
                            good = False
                            break
                    if not good:
                        break
                if not good:
                    break

            if good:
                if (perm1, udiff1, mul_val1) not in ret_list:
                    ret_list[(perm1, udiff1, mul_val1)] = []
                ret_list[(perm1, udiff1, mul_val1)] += [(perm2, udiff2, mul_val2)]
    return ret_list


def will_formula_work(u, v):
    u, v = Permutation(u), Permutation(v)
    muv = uncode(theta(v))
    vn1muv = (~v) * muv
    while inv(vn1muv) > 0:
        found_one = False
        for i in range(len(vn1muv) - 1):
            if vn1muv[i] > vn1muv[i + 1]:
                found_one = True
                if i < len(u) - 1 and u[i] > u[i + 1]:
                    return False
                # vn1muv[i], vn1muv[i + 1] = vn1muv[i + 1], vn1muv[i]
                vn1muv = vn1muv.swap(i, i + 1)
                break
        if found_one:
            return False
    return True


def try_reduce_u(u, v, w):
    if one_dominates(u, w):
        return u, v, w
    u2 = u
    v2 = v
    w2 = w
    cu = code(u)
    for i in range(len(u2) - 2, -1, -1):
        if cu[i] == 0 and i < len(cu) - 1 and cu[i + 1] != 0:
            if i >= len(v2) - 1 or v2[i] < v2[i + 1]:
                # u2[i], u2[i + 1] = u2[i + 1], u2[i]
                # if i > len(w2) - 1:
                #     w2 += list(range(len(w2) + 1, i + 3))
                # w2[i + 1], w2[i] = w2[i], w2[i + 1]
                u2 = u2.swap(i, i + 1)
                w2 = w2.swap(i, i + 1)
                if one_dominates(u2, w):
                    # return Permutation(u2), Permutation(v2), Permutation(w2)
                    return u2, v2, w2
                return try_reduce_u(u2, v2, w2)
            if i < len(w2) - 1 and w2[i] > w2[i + 1]:
                # ERROR?
                # u2[i], u2[i + 1] = u2[i + 1], u2[i]
                # v2[i], v2[i + 1] = v2[i + 1], v2[i]
                u2 = u2.swap(i, i + 1)
                v2 = v2.swap(i, i + 1)
                return try_reduce_u(u2, v2, w2)
            # return Permutation(u2), Permutation(v2), Permutation(w2)
            return u2, v2, w2
    # return Permutation(u2), Permutation(v2), Permutation(w2)
    return u2, v2, w2


def reduce_descents(u, v, w):
    found_one = True
    u2 = Permutation(u)
    v2 = Permutation(v)
    w2 = Permutation(w)
    while found_one:
        found_one = False
        if will_formula_work(u2, v2) or will_formula_work(v2, u2) or one_dominates(u2, w2) or is_reducible(v2) or inv(w2) - inv(u2) == 1:
            break
        for i in range(len(w2) - 2, -1, -1):
            if w2[i] > w2[i + 1] and i < len(v2) - 1 and v2[i] > v2[i + 1] and (i >= len(u2) - 1 or u2[i] < u2[i + 1]):
                w2 = w2.swap(i, i + 1)
                v2 = v2.swap(i, i + 1)
                found_one = True
            elif w2[i] > w2[i + 1] and i < len(u2) - 1 and u2[i] > u2[i + 1] and (i >= len(v2) - 1 or v2[i] < v2[i + 1]):
                w2 = w2.swap(i, i + 1)
                u2 = u2.swap(i, i + 1)
                found_one = True
            if found_one:
                break
    # return Permutation(u2), Permutation(v2), Permutation(w2)
    return u2, v2, w2


def is_reducible(v):
    c03 = code(v)
    found0 = False
    good = True
    for i in range(len(c03)):
        if c03[i] == 0:
            found0 = True
        elif c03[i] != 0 and found0:
            good = False
            break
    return good


#
def try_reduce_v(u, v, w):
    if is_reducible(v):
        return u, v, w
    u2 = u
    v2 = v
    w2 = w
    cv = code(v2)
    for i in range(len(v2) - 2, -1, -1):
        if cv[i] == 0 and i < len(cv) - 1 and cv[i + 1] != 0:
            if i >= len(u2) - 1 or u2[i] < u2[i + 1]:
                # v2[i], v2[i + 1] = v2[i + 1], v2[i]
                v2 = v2.swap(i, i + 1)
                # if i >= len(w2) - 1:
                #     w2 += list(range(len(w2) + 1, i + 3))
                # w2[i + 1], w2[i] = w2[i], w2[i + 1]
                w2 = w2.swap(i, i + 1)
                if is_reducible(v2):
                    return Permutation(u2), Permutation(v2), Permutation(w2)
                return try_reduce_v(u2, v2, w2)
            if i < len(w2) - 1 and w2[i] > w2[i + 1]:
                # u2[i], u2[i + 1] = u2[i + 1], u2[i]
                # v2[i], v2[i + 1] = v2[i + 1], v2[i]
                u2 = u2.swap(i, i + 1)
                v2 = v2.swap(i, i + 1)
                return try_reduce_v(u2, v2, w2)
            # return Permutation(u2), Permutation(v2), Permutation(w2)
            return u2, v2, w2
    return u2, v2, w2


def reduce_coeff(u, v, w):
    t_mu_u_t = theta(inverse(u))
    t_mu_v_t = theta(inverse(v))

    mu_u_inv = uncode(t_mu_u_t)
    mu_v_inv = uncode(t_mu_v_t)

    t_mu_u = p_trans(t_mu_u_t)
    t_mu_v = p_trans(t_mu_v_t)

    t_mu_u += [0 for i in range(len(t_mu_u), max(len(t_mu_u), len(t_mu_v)))]
    t_mu_v += [0 for i in range(len(t_mu_v), max(len(t_mu_u), len(t_mu_v)))]

    t_mu_uv = [t_mu_u[i] + t_mu_v[i] for i in range(len(t_mu_u))]
    t_mu_uv_t = p_trans(t_mu_uv)

    mu_uv_inv = uncode(t_mu_uv_t)

    if inv(w * mu_uv_inv) != inv(mu_uv_inv) - inv(w):
        return u, v, w

    # umu = mulperm(list(u), mu_u_inv)
    # vmu = mulperm(list(v), mu_v_inv)
    # wmu = mulperm(list(w), mu_uv_inv)
    umu = u * mu_u_inv
    vmu = v * mu_v_inv
    wmu = w * mu_uv_inv

    t_mu_w = theta(inverse(wmu))

    mu_w = uncode(t_mu_w)

    w_prime = wmu * mu_w

    if w == w_prime:
        return u, v, w
    # if permtrim(list(w)) == permtrim(w_prime):
    #     return (permtrim(list(u)), permtrim(list(v)), permtrim(list(w)))

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

    mu_w_A = uncode(mu_A(code(mu_w), A))
    mu_w_B = uncode(mu_A(code(mu_w), B))

    return (umu * mu_w_A, vmu * mu_w_B, w_prime)


def pull_out_var(vnum, v):
    v = Permutation(v)
    vup = v
    if vnum >= len(v):
        return {((), v)}
    vpm_list = [(vup, 0)]
    ret_list = set()
    for p in range(len(v) + 1 - vnum):
        vpm_list2 = []
        for vpm, b in vpm_list:
            if vpm[vnum - 1] == len(v) + 1:
                vpm2 = [*vpm]
                vpm2.pop(vnum - 1)
                vp = permtrim(vpm2)
                ret_list.add(
                    (
                        tuple([v[i] for i in range(vnum, len(v)) if ((i > len(vp) and v[i] == i) or (i <= len(vp) and v[i] == vp[i - 1]))]),
                        vp,
                    ),
                )
            for j in range(vnum, len(vup) + 2):
                if vpm[j] <= b:
                    continue
                for i in range(vnum):
                    if has_bruhat_ascent(vpm, i, j):
                        vpm_list2 += [(vpm.swap(i, j), vpm[j])]
        vpm_list = vpm_list2
    for vpm, b in vpm_list:
        if vpm[vnum - 1] == len(v) + 1:
            vpm2 = [*vpm]
            vpm2.pop(vnum - 1)
            vp = permtrim(vpm2)
            ret_list.add(
                (
                    tuple([v[i] for i in range(vnum, len(v)) if ((i > len(vp) and v[i] == i) or (i <= len(vp) and v[i] == vp[i - 1]))]),
                    vp,
                ),
            )

    return ret_list


def divdiffable(v, u):
    inv_v = inv(v)
    inv_u = inv(u)
    perm2 = v * (~u)
    if inv(perm2) != inv_v - inv_u:
        return []
    return perm2


# perms and inversion diff
# def kdown_perms(perm, monoperm, p, k):
#     perm = Permutation(perm)
#     monoperm = Permutation(monoperm)
#     imonoperm = ~monoperm
#     inv_m = inv(monoperm)
#     inv_p = inv(perm)
#     full_perm_list = []
#     # perm = Permutation(perm)
#     # print(f"{perm=} {monoperm=} {inv_m=} {inv_p=} {perm*monoperm=}")
#     if inv(perm * monoperm) == inv_m - inv_p:
#         full_perm_list += [(perm, 0, 1)]
#     down_perm_list = [(perm, 1, perm * monoperm)]
#     if len(perm) < k:
#         return full_perm_list
#     a2 = k - 1
#     for pp in range(1, p + 1):
#         down_perm_list2 = []
#         for perm2, s, test_perm in down_perm_list:
#             L = len(perm2)
#             inv_test_perm = inv(test_perm)
#             if k > L:
#                 continue
#             s2 = -s
#             for b in chain(range(k - 1), range(k, L)):
#                 if perm2[b] != perm[b]:
#                     continue
#                 if b < a2:
#                     i, j = b, a2
#                 else:
#                     i, j, s2 = a2, b, s
#                 # print(f"{perm2=} {i=} {j=}")
#                 if has_bruhat_descent(perm2, i, j):
#                     new_perm = perm2.swap(i, j)
#                     down_perm_list2 += [(new_perm, s2, test_perm.swap(imonoperm[i] - 1,imonoperm[j] - 1))]
#                     add_inv = 1 if test_perm[imonoperm[i] - 1] < test_perm[imonoperm[j] - 1] else -1
#                     if inv_test_perm + add_inv == inv_m - inv_p + pp:
#                         full_perm_list += [(new_perm, pp, s2)]
#         down_perm_list = down_perm_list2
#     return full_perm_list


def kdown_perms(perm, monoperm, p, k):
    perm = Permutation(perm)
    monoperm = Permutation(monoperm)
    inv_m = inv(monoperm)
    inv_p = inv(perm)
    full_perm_list = []
    if inv(perm * monoperm) == inv_m - inv_p:
        full_perm_list += [(perm, 0, 1)]

    down_perm_list = [(perm, S.One)]
    if len(perm) < k:
        return full_perm_list
    a2 = k - 1
    for pp in range(1, p + 1):
        down_perm_list2 = []
        g_inv = inv_m - inv_p + pp
        for perm2, s in down_perm_list:
            L = len(perm2)
            if k > L:
                continue
            s2 = -s
            rg = [i for i in range(k - 1) if perm2[i] == perm[i]]
            for b in rg:
                if has_bruhat_descent(perm2, b, a2):
                    new_perm = perm2.swap(b, a2)
                    down_perm_list2 += [(new_perm, s2)]
                    if inv(new_perm * monoperm) == g_inv:
                        full_perm_list += [(new_perm, pp, s2)]
            rg = [i for i in range(k, L) if perm2[i] == perm[i]]
            for b in rg:
                if has_bruhat_descent(perm2, a2, b):
                    new_perm = perm2.swap(a2, b)
                    down_perm_list2 += [(new_perm, s)]
                    if inv(new_perm * monoperm) == g_inv:
                        full_perm_list += [(new_perm, pp, s)]
        down_perm_list = down_perm_list2
    return full_perm_list

def rc_graph_set(perm):
    if perm.inv == 0:
        return {((),())}
    ret = set()
    L = pull_out_var(1, perm)
    for index_list, new_perm in L:
        rc_set = rc_graph_set(new_perm)
        lsort = sorted(index_list, reverse=True)
        for labels, word in rc_set:
            new_labels = tuple(([1]*len(index_list)) + [label + 1 for label in labels])
            new_word = tuple(lsort+[word_s + 1 for word_s in word])
            ret.add((new_labels, new_word))
    return ret


@cache
def compute_vpathdicts_cached(th, vmu):
    vpathdicts = [{} for index in range(len(th))]
    vpathdicts[-1][vmu] = None
    thL = len(th)

    top = code(~Permutation(uncode(th)))
    for i in range(thL - 1, -1, -1):
        top2 = code(~Permutation(uncode(top)))
        while top2[-1] == 0:
            top2.pop()
        top2.pop()
        top = code(~Permutation(uncode(top2)))
        # print(f"{top=}")
        monoperm = Permutation(uncode(top))
        # print(f"{monoperm=}")
        k = i + 1
        for last_perm in vpathdicts[i]:
            newperms = kdown_perms(last_perm, monoperm, th[i], k)
            vpathdicts[i][last_perm] = newperms
            if i > 0:
                for trip in newperms:
                    vpathdicts[i - 1][trip[0]] = None
    vpathdicts2 = [{} for i in range(len(th))]
    for i in range(len(th)):
        for key, valueset in vpathdicts[i].items():
            for value in valueset:
                key2 = value[0]
                if key2 not in vpathdicts2[i]:
                    vpathdicts2[i][key2] = set()
                v2 = value[2]
                vpathdicts2[i][key2].add((key, value[1], v2))
    return vpathdicts2


def compute_vpathdicts(th, vmu):
    return compute_vpathdicts_cached(tuple(th), vmu)


def check_blocks(qv, parabolic_index):
    blocks = []
    cur_block = []
    last_val = -1
    for i in range(len(parabolic_index)):
        if last_val + 1 != parabolic_index[i] and i != 0:
            blocks += [cur_block]
            cur_block = []
        last_val = parabolic_index[i]
        cur_block += [last_val]
    if len(cur_block) != 0:
        blocks += [cur_block]
    # print(f"{parabolic_index=}")
    # print(f"{blocks=}")
    for block in blocks:
        for i in range(len(block)):
            for j in range(i, len(block)):
                val = 0
                for k in range(i, j + 1):
                    val += omega(block[k], qv)
                if val != 0 and val != -1:
                    return False
    return True


def reduce_q_coeff(u, v, w, qv):
    for i in range(len(qv)):
        if sg(i, v) == 1 and sg(i, u) == 0 and sg(i, w) + omega(i + 1, qv) == 1:
            ret_v = v.swap(i, i + 1)
            ret_w = w.swap(i, i + 1)
            qv_ret = [*qv]
            if sg(i, w) == 0:
                qv_ret[i] -= 1
            return u, ret_v, ret_w, qv_ret, True
        if (sg(i, u) == 1 and sg(i, v) == 0 and sg(i, w) + omega(i + 1, qv) == 1) or (sg(i, u) == 1 and sg(i, v) == 1 and sg(i, w) + omega(i + 1, qv) == 2):
            ret_u = u.swap(i, i + 1)
            ret_w = w.swap(i, i + 1)
            qv_ret = [*qv]
            if sg(i, w) == 0:
                qv_ret[i] -= 1
            return ret_u, v, ret_w, qv_ret, True
    return u, v, w, qv, False


def reduce_q_coeff_u_only(u, v, w, qv):
    for i in range(len(qv)):
        if (sg(i, u) == 1 and sg(i, v) == 0 and sg(i, w) + omega(i + 1, qv) == 1) or (sg(i, u) == 1 and sg(i, v) == 1 and sg(i, w) + omega(i + 1, qv) == 2):
            # ret_u = [*u]
            # ret_u[i], ret_u[i + 1] = ret_u[i + 1], ret_u[i]
            ret_u = u.swap(i, i + 1)
            # ret_w = [*w] + [j + 1 for j in range(len(w), i + 2)]
            # ret_w[i], ret_w[i + 1] = ret_w[i + 1], ret_w[i]
            ret_w = w.swap(i, i + 1)
            qv_ret = [*qv]
            if sg(i, w) == 0:
                qv_ret[i] -= 1
            return ret_u, v, ret_w, qv_ret, True
    return u, v, w, qv, False


def elem_sym_perms_q(orig_perm, p, k, q_var=q_var):
    total_list = [(orig_perm, 0, 1)]
    up_perm_list = [(orig_perm, 1, 1000)]
    for pp in range(p):
        perm_list = []
        for up_perm, val, last_j in up_perm_list:
            pos_list = [i for i in range(k) if up_perm[i] == orig_perm[i]]
            for j in range(min(max(k + 1, len(up_perm)), last_j), k - 1, -1):
                for i in pos_list:
                    ct = count_bruhat(up_perm, i, j)
                    if ct == 1 or ct == 2 * (i - j) + 1:
                        new_perm_add = up_perm.swap(i, j)
                        # print(f"{new_perm_add=}")
                        new_val = val
                        if ct < 0:
                            new_val *= prod([q_var[index] for index in range(i + 1, j + 1)])
                        perm_list += [(new_perm_add, new_val, j)]
                        total_list += [(new_perm_add, pp + 1, new_val)]
        up_perm_list = perm_list
    return total_list


def elem_sym_perms_q_op(orig_perm, p, k, n, q_var=q_var):
    total_list = [(orig_perm, 0, 1)]
    up_perm_list = [(orig_perm, 1, k)]
    for pp in range(p):
        perm_list = []
        for up_perm, val, last_j in up_perm_list:
            up_perm2 = [*up_perm]
            if len(up_perm) < n:
                up_perm2 += [i + 1 for i in range(len(up_perm2), n)]
            pos_list = [i for i in range(k) if (i >= len(orig_perm) and up_perm2[i] == i + 1) or (i < len(orig_perm) and up_perm2[i] == orig_perm[i])]
            for j in range(last_j, n):
                for i in pos_list:
                    ct = count_bruhat(up_perm2, i, j)
                    # print(f"{up_perm2=} {ct=} {i=} {j=} {k=} {pp=}")
                    if ct == -1 or ct == 2 * (j - i) - 1:
                        new_perm = [*up_perm2]
                        new_perm[i], new_perm[j] = new_perm[j], new_perm[i]
                        new_perm_add = Permutation(new_perm)
                        new_val = val
                        if ct > 0:
                            new_val *= prod([q_var[index] for index in range(i + 1, j + 1)])
                        perm_list += [(new_perm_add, new_val, j)]
                        total_list += [(new_perm_add, pp + 1, new_val)]
        up_perm_list = perm_list
    return total_list


def elem_sym_perms(orig_perm, p, k):
    orig_perm = Permutation(orig_perm)
    total_list = [(orig_perm, 0)]
    up_perm_list = [(orig_perm, 1000000000)]
    for pp in range(p):
        perm_list = []
        for up_perm, last in up_perm_list:
            pos_list = [i for i in range(k) if up_perm[i] < last]
            for j in range(k, max(k + 2, len(up_perm) + 1)):
                if up_perm[j] >= last:
                    continue
                for i in pos_list:
                    if has_bruhat_ascent(up_perm, i, j):
                        new_perm_add = up_perm.swap(i, j)
                        perm_list += [(new_perm_add, up_perm[j])]
                        total_list += [(new_perm_add, pp + 1)]
        up_perm_list = perm_list
    return total_list


def elem_sym_positional_perms(orig_perm, p, *k):
    k = {i - 1 for i in k}
    orig_perm = Permutation(orig_perm)
    total_list = {(orig_perm, 0, 1)}
    up_perm_list = {(orig_perm, 1)}
    for pp in range(p):
        perm_list = set()
        for up_perm, sign in up_perm_list:
            # pos_list = [i for i in range(k) if up_perm[i] < last]
            rg = [q for q in range(len(up_perm) + max(k) + 1) if q not in k]
            pos_list = [i for i in k if up_perm[i] == orig_perm[i]]
            for i in pos_list:
                for j in rg:
                    a, b = (i, j) if i < j else (j, i)
                    if has_bruhat_ascent(up_perm, a, b):
                        new_perm_add = up_perm.swap(a, b)
                        new_sign = sign if i < j else -sign
                        perm_list.add((new_perm_add, new_sign))
                        total_list.add((new_perm_add, pp + 1, new_sign))
        up_perm_list = perm_list
    return total_list


# def elem_sym_perms_q(orig_perm, p, k, q_var=q_var):
#     total_list = [(orig_perm, 0, 1)]
#     up_perm_list = [(orig_perm, 1, 1000)]
#     for pp in range(p):
#         perm_list = []
#         for up_perm, val, last_j in up_perm_list:
#             pos_list = [i for i in range(k) if up_perm[i] == orig_perm[i]]
#             for j in range(min(max(k + 1, len(up_perm)), last_j), k - 1, -1):
#                 for i in pos_list:
#                     ct = count_bruhat(up_perm, i, j)
#                     if ct == 1 or ct == 2 * (i - j) + 1:
#                         new_perm_add = up_perm.swap(i, j)
#                         # print(f"{new_perm_add=}")
#                         new_val = val
#                         if ct < 0:
#                             new_val *= prod([q_var[index] for index in range(i + 1, j + 1)])
#                         perm_list += [(new_perm_add, new_val, j)]
#                         total_list += [(new_perm_add, pp + 1, new_val)]
#         up_perm_list = perm_list
#     return total_list


# def elem_sym_positional_perms_q(orig_perm, p, *k, q_var=q_var):
#     k = {i - 1 for i in k}
#     orig_perm = Permutation(orig_perm)
#     total_list = {(orig_perm, 0, S.One, S.One)}
#     up_perm_list = {(orig_perm, S.One, S.One, 1000)}
#     for pp in range(p):
#         perm_list = set()
#         for up_perm, sign, val, last_j in up_perm_list:
#             pos_list = [i for i in k if up_perm[i] == orig_perm[i]]
#             rg = [q for q in range(min(len(up_perm) + max(k), last_j + 1)) if q not in k]
#             for j in rg:
#                 for i in pos_list:
#                     a, b = (i, j) if i < j else (j, i)
#                     ct = count_bruhat(up_perm, a, b)
#                     if ct == 1 or ct == 2 * (a - b) + 1:
#                         new_perm_add = up_perm.swap(i, j)
#                         # print(f"{new_perm_add=} {i=} {j=} {p=} {last_j=} {pos_list=}")
#                         new_sign = sign if i < j else -sign
#                         new_val = val
#                         if ct < 0:
#                             new_val *= prod([q_var[index] for index in range(a + 1, b + 1)])
#                         # print(f"{val=} {new_val=}")
#                         perm_list.add((new_perm_add, new_sign, new_val, j))
#                         total_list.add((new_perm_add, pp + 1, new_sign, new_val))
#         up_perm_list = perm_list
#     return total_list


# this is a quantization
def elem_sym_positional_perms_q(orig_perm, p, *k, q_var=q_var):
    k = {i - 1 for i in k}
    # print(f"{k=}")
    # sorting perm
    # lookup_mask = Permutation.sorting_perm(k)
    orig_perm = Permutation(orig_perm)
    total_list = {(orig_perm, 0, S.One, S.One)}
    up_perm_list = {(orig_perm, S.One, S.One, 1000)}
    for pp in range(p):
        perm_list = set()
        for up_perm, sign, val, last_j in up_perm_list:
            # print(f"{up_perm=}")
            pos_list = [i for i in k if up_perm[i] == orig_perm[i]]
            rg = [q for q in range(min(last_j, len(up_perm) + max(k)), -1, -1) if q not in k]
            # print(f"{rg=}")
            # print(f"{pos_list=}")
            for i in pos_list:
                for j in rg:
                    a, b = (i, j) if i < j else (j, i)
                    ct = count_bruhat(up_perm, a, b)
                    # print(f"{up_perm=} {a=} {b=} {ct=}")
                    if ct == 1 or ct == 2 * (a - b) + 1:
                        new_perm_add = up_perm.swap(i, j)
                        # print(f"{new_perm_add=} {i=} {j=} {p=} {pos_list=}")
                        new_sign = sign if i < j else -sign
                        new_val = val
                        if ct < 0:
                            new_val *= prod([q_var[index] for index in range(a + 1, b + 1)])
                        # print(f"{val=} {new_val=}")
                        perm_list.add((new_perm_add, new_sign, new_val, j))
                        total_list.add((new_perm_add, pp + 1, new_sign, new_val))
        up_perm_list = perm_list
    return total_list


def complete_sym_perms_op(orig_perm, p, k):
    from schubmult.utils.perm_utils import has_bruhat_descent

    orig_perm = Permutation(orig_perm)
    total_dict = {orig_perm: []}
    up_perm_list = [(orig_perm, 1000000000)]
    for pp in range(p):
        perm_list = []
        new_total_dict = {}
        for up_perm, last in up_perm_list:
            pos_list = [j for j in range(len(orig_perm)) if up_perm[j] == orig_perm[j]]
            for i in range(k):
                for j in pos_list:
                    if has_bruhat_descent(up_perm, i, j):
                        new_perm_add = up_perm.swap(i, j)
                        perm_list += [(new_perm_add, up_perm[j])]
                        new_total_dict[new_perm_add] = [(i + 1, j + 1)] + total_dict[up_perm]
        up_perm_list = perm_list
        total_dict = new_total_dict
    return total_dict

def complete_sym_perms(orig_perm, p, k):
    from schubmult.utils.perm_utils import has_bruhat_ascent

    orig_perm = Permutation(orig_perm)
    total_dict = {orig_perm: []}
    up_perm_list = [(orig_perm, 1000000000)]
    for pp in range(p):
        perm_list = []
        new_total_dict = {}
        for up_perm, last in up_perm_list:
            pos_list = [j for j in range(len(orig_perm) + 10) if up_perm[j] == orig_perm[j]]
            for i in range(k):
                for j in pos_list:
                    if has_bruhat_ascent(up_perm, i, j):
                        new_perm_add = up_perm.swap(i, j)
                        perm_list += [(new_perm_add, up_perm[j])]
                        new_total_dict[new_perm_add] = total_dict[up_perm] + [(i + 1, j + 1)]
        up_perm_list = perm_list
        total_dict = new_total_dict
    return total_dict


def complete_sym_positional_perms(orig_perm, p, *k):
    k = {i - 1 for i in k}
    orig_perm = Permutation(orig_perm)
    total_list = {(orig_perm, 0, 1)}
    up_perm_list = {(orig_perm, 1)}
    # print(f"{orig_perm=}")
    for pp in range(p):
        perm_list = set()
        for up_perm, sign in up_perm_list:
            # pos_list = [i for i in range(k) if up_perm[i] < last]
            rg = [q for q in range(len(up_perm) + max(k) + 1) if q not in k and up_perm[q] == orig_perm[q]]
            for j in rg:
                for i in k:
                    a, b = (i, j) if i < j else (j, i)
                    if has_bruhat_ascent(up_perm, a, b):
                        # print(f"bruhat ascent on {up_perm=} at {(a,b)=}")
                        new_perm_add = up_perm.swap(a, b)
                        # print(f"{up_perm.inv - orig_perm.inv=}")
                        # print(f"{new_perm_add.inv - orig_perm.inv=}")
                        new_sign = sign if i < j else -sign
                        perm_list.add((new_perm_add, new_sign))
                        total_list.add((new_perm_add, pp + 1, new_sign))
        up_perm_list = perm_list
    return total_list

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
    return total_list


def elem_sym_perms_op(orig_perm, p, k):
    total_list = [(orig_perm, 0)]
    up_perm_list = [(orig_perm, k)]
    for pp in range(p):
        perm_list = []
        for up_perm, last in up_perm_list:
            up_perm2 = [*up_perm]
            if len(up_perm2) < k + 1:
                up_perm2 += [i + 1 for i in range(len(up_perm2), k + 2)]
            pos_list = [i for i in range(k) if getpermval(up_perm2, i) == getpermval(orig_perm, i)]
            for j in range(last, len(up_perm2)):
                for i in pos_list:
                    if has_bruhat_descent(up_perm2, i, j):
                        new_perm_add = up_perm.swap(i, j)
                        perm_list += [(new_perm_add, j)]
                        total_list += [(new_perm_add, pp + 1)]
        up_perm_list = perm_list
    return total_list


def is_split_two(u, v, w):
    if inv(w) - inv(u) != 2:
        return False, []
    diff_perm = (~v) * w
    identity = [i + 1 for i in range(len(diff_perm))]
    cycles = []
    for i in range(len(identity)):
        if diff_perm[i] != identity[i]:
            cycle0 = set()
            cycle = {i + 1}
            last = i
            while len(cycle0) != len(cycle):
                cycle0 = cycle
                last = diff_perm[last] - 1
                cycle.add(last + 1)
            if len(cycle) > 1 and cycle not in cycles:
                cycles += [cycle]
            if len(cycles) > 2:
                break
    if len(cycles) == 2:
        return True, cycles
    return False, []


def is_coeff_irreducible(u, v, w):
    return (
        not will_formula_work(u, v)
        and not will_formula_work(v, u)
        and not one_dominates(u, w)
        and not is_reducible(v)
        and inv(w) - inv(u) > 1
        and not is_split_two(u, v, w)[0]
        and len([i for i in code(v) if i != 0]) > 1
    )


def is_hook(cd):
    started = False
    done = False
    found_zero_after = False
    for i in range(len(cd)):
        if (done or found_zero_after) and cd[i] != 0:
            return False
        if cd[i] == 1 and not started:
            started = True
        if cd[i] > 1:
            done = True
        if started and cd[i] == 0:
            found_zero_after = True
    if started or done:
        return True
    return False
