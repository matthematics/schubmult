from bisect import bisect_left
from functools import cache
from itertools import chain

import numpy as np
from symengine import Mul, Pow, symarray, sympify

zero = sympify(0)
n = 100

q_var = symarray("q", n)


def getpermval(perm, index):
    if index < len(perm):
        return perm[index]
    return index + 1


def inv(perm):
    L = len(perm)
    v = list(range(1, L + 1))
    ans = 0
    for i in range(L):
        itr = bisect_left(v, perm[i])
        ans += itr
        v = v[:itr] + v[itr + 1 :]
    return ans


def code(perm):
    L = len(perm)
    ret = []
    v = list(range(1, L + 1))
    for i in range(L - 1):
        itr = bisect_left(v, perm[i])
        ret += [itr]
        v = v[:itr] + v[itr + 1 :]
    return ret


def mulperm(perm1, perm2):
    if len(perm1) < len(perm2):
        return [perm1[perm2[i] - 1] if perm2[i] <= len(perm1) else perm2[i] for i in range(len(perm2))]
    return [perm1[perm2[i] - 1] for i in range(len(perm2))] + perm1[len(perm2) :]


def uncode(cd):
    cd2 = [*cd]
    if cd2 == []:
        return [1, 2]
    max_required = max([cd2[i] + i for i in range(len(cd2))])
    cd2 += [0 for i in range(len(cd2), max_required)]
    fullperm = [i + 1 for i in range(len(cd2) + 1)]
    perm = []
    for i in range(len(cd2)):
        perm += [fullperm.pop(cd2[i])]
    perm += [fullperm[0]]
    return perm


def reversecode(perm):
    ret = []
    for i in range(len(perm) - 1, 0, -1):
        ret = [0, *ret]
        for j in range(i, -1, -1):
            if perm[i] > perm[j]:
                ret[-1] += 1
    return ret


def reverseuncode(cd):
    cd2 = list(cd)
    if cd2 == []:
        return [1, 2]
    # max_required = max([cd2[i]+i for i in range(len(cd2))])
    # cd2 += [0 for i in range(len(cd2),max_required)]
    fullperm = [i + 1 for i in range(len(cd2) + 1)]
    perm = []
    for i in range(len(cd2) - 1, 0, -1):
        perm = [fullperm[cd2[i]], *perm]
        fullperm.pop(cd2[i])
    perm += [fullperm[0]]
    return perm


def inverse(perm):
    retperm = [0 for i in range(len(perm))]
    for i in range(len(perm)):
        retperm[perm[i] - 1] = i + 1
    return retperm


def permtrim(perm):
    L = len(perm)
    while L > 2 and perm[-1] == L:
        L = perm.pop() - 1
    return perm


def has_bruhat_descent(perm, i, j):
    if perm[i] < perm[j]:
        return False
    for p in range(i + 1, j):
        if perm[i] > perm[p] and perm[p] > perm[j]:
            return False
    return True


def count_bruhat(perm, i, j):
    up_amount = 0
    if perm[i] < perm[j]:
        up_amount = 1
    else:
        up_amount = -1
    for k in range(i + 1, j):
        if perm[i] < perm[k] and perm[k] < perm[j]:
            up_amount += 2
        elif perm[i] > perm[k] and perm[k] > perm[j]:
            up_amount -= 2
    return up_amount


def has_bruhat_ascent(perm, i, j):
    if perm[i] > perm[j]:
        return False
    for p in range(i + 1, j):
        if perm[i] < perm[p] and perm[p] < perm[j]:
            return False
    return True


def elem_sym_perms(orig_perm, p, k):
    total_list = [(orig_perm, 0)]
    up_perm_list = [(orig_perm, 1000000000)]
    for pp in range(p):
        perm_list = []
        for up_perm, last in up_perm_list:
            up_perm2 = [*up_perm, len(up_perm) + 1]
            if len(up_perm2) < k + 1:
                up_perm2 += [i + 1 for i in range(len(up_perm2), k + 2)]
            pos_list = [i for i in range(k) if up_perm2[i] < last]
            for j in range(k, len(up_perm2)):
                if up_perm2[j] >= last:
                    continue
                for i in pos_list:
                    if has_bruhat_ascent(up_perm2, i, j):
                        new_perm = [*up_perm2]
                        new_perm[i], new_perm[j] = new_perm[j], new_perm[i]
                        if new_perm[-1] == len(new_perm):
                            new_perm_add = tuple(new_perm[:-1])
                        else:
                            new_perm_add = tuple(new_perm)
                        perm_list += [(new_perm_add, up_perm2[j])]
                        total_list += [(new_perm_add, pp + 1)]
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
                        new_perm = [*up_perm2]
                        new_perm[i], new_perm[j] = new_perm[j], new_perm[i]
                        new_perm_add = tuple(permtrim(new_perm))
                        perm_list += [(new_perm_add, j)]
                        total_list += [(new_perm_add, pp + 1)]
        up_perm_list = perm_list
    return total_list


def strict_theta(u):
    ret = [*trimcode(u)]
    did_one = True
    while did_one:
        did_one = False
        for i in range(len(ret) - 2, -1, -1):
            if ret[i + 1] != 0 and ret[i] <= ret[i + 1]:
                ret[i], ret[i + 1] = ret[i + 1] + 1, ret[i]
                did_one = True
                break
    while len(ret) > 0 and ret[-1] == 0:
        ret.pop()
    return ret


def elem_sym_perms_q(orig_perm, p, k, q_var=q_var):
    total_list = [(orig_perm, 0, 1)]
    up_perm_list = [(orig_perm, 1, 1000)]
    for pp in range(p):
        perm_list = []
        for up_perm, val, last_j in up_perm_list:
            up_perm2 = [*up_perm, len(up_perm) + 1]
            if len(up_perm2) < k + 1:
                up_perm2 += [i + 1 for i in range(len(up_perm2), k + 2)]
            pos_list = [i for i in range(k) if (i >= len(orig_perm) and up_perm2[i] == i + 1) or (i < len(orig_perm) and up_perm2[i] == orig_perm[i])]
            for j in range(min(len(up_perm2) - 1, last_j), k - 1, -1):
                for i in pos_list:
                    ct = count_bruhat(up_perm2, i, j)
                    # print(f"{up_perm2=} {ct=} {i=} {j=} {k=} {pp=}")
                    if ct == 1 or ct == 2 * (i - j) + 1:
                        new_perm = [*up_perm2]
                        new_perm[i], new_perm[j] = new_perm[j], new_perm[i]
                        new_perm_add = tuple(permtrim(new_perm))
                        new_val = val
                        if ct < 0:
                            new_val *= np.prod([q_var[index] for index in range(i + 1, j + 1)])
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
                        new_perm_add = tuple(permtrim(new_perm))
                        new_val = val
                        if ct > 0:
                            new_val *= np.prod([q_var[index] for index in range(i + 1, j + 1)])
                        perm_list += [(new_perm_add, new_val, j)]
                        total_list += [(new_perm_add, pp + 1, new_val)]
        up_perm_list = perm_list
    return total_list


def q_vector(q_exp, q_var=q_var):
    qvar_list = q_var.tolist()
    ret = []

    if q_exp == 1:
        return ret
    if q_exp in q_var:
        i = qvar_list.index(q_exp)
        return [0 for j in range(i - 1)] + [1]
    if isinstance(q_exp, Pow):
        qv = q_exp.args[0]
        expon = int(q_exp.args[1])
        i = qvar_list.index(qv)
        return [0 for j in range(i - 1)] + [expon]
    if isinstance(q_exp, Mul):
        for a in q_exp.args:
            v1 = q_vector(a)
            v1 += [0 for i in range(len(v1), len(ret))]
            ret += [0 for i in range(len(ret), len(v1))]
            ret = [ret[i] + v1[i] for i in range(len(ret))]
        return ret

    return None


def omega(i, qv):
    i = i - 1
    if len(qv) == 0 or i > len(qv):
        return 0
    if i == 0:
        if len(qv) == 1:
            return 2 * qv[0]
        return 2 * qv[0] - qv[1]
    if i == len(qv):
        return -qv[-1]
    if i == len(qv) - 1:
        return 2 * qv[-1] - qv[-2]
    return 2 * qv[i] - qv[i - 1] - qv[i + 1]


def sg(i, w):
    if i >= len(w) - 1 or w[i] < w[i + 1]:
        return 0
    return 1


def reduce_q_coeff(u, v, w, qv):
    for i in range(len(qv)):
        if sg(i, v) == 1 and sg(i, u) == 0 and sg(i, w) + omega(i + 1, qv) == 1:
            ret_v = [*v]
            ret_v[i], ret_v[i + 1] = ret_v[i + 1], ret_v[i]
            ret_w = [*w] + [j + 1 for j in range(len(w), i + 2)]
            ret_w[i], ret_w[i + 1] = ret_w[i + 1], ret_w[i]
            qv_ret = [*qv]
            if sg(i, w) == 0:
                qv_ret[i] -= 1
            return u, tuple(permtrim(ret_v)), tuple(permtrim(ret_w)), qv_ret, True
        if (sg(i, u) == 1 and sg(i, v) == 0 and sg(i, w) + omega(i + 1, qv) == 1) or (sg(i, u) == 1 and sg(i, v) == 1 and sg(i, w) + omega(i + 1, qv) == 2):
            ret_u = [*u]
            ret_u[i], ret_u[i + 1] = ret_u[i + 1], ret_u[i]
            ret_w = [*w] + [j + 1 for j in range(len(w), i + 2)]
            ret_w[i], ret_w[i + 1] = ret_w[i + 1], ret_w[i]
            qv_ret = [*qv]
            if sg(i, w) == 0:
                qv_ret[i] -= 1
            return tuple(permtrim(ret_u)), v, tuple(permtrim(ret_w)), qv_ret, True
    return u, v, w, qv, False


def reduce_q_coeff_u_only(u, v, w, qv):
    for i in range(len(qv)):
        if (sg(i, u) == 1 and sg(i, v) == 0 and sg(i, w) + omega(i + 1, qv) == 1) or (sg(i, u) == 1 and sg(i, v) == 1 and sg(i, w) + omega(i + 1, qv) == 2):
            ret_u = [*u]
            ret_u[i], ret_u[i + 1] = ret_u[i + 1], ret_u[i]
            ret_w = [*w] + [j + 1 for j in range(len(w), i + 2)]
            ret_w[i], ret_w[i + 1] = ret_w[i + 1], ret_w[i]
            qv_ret = [*qv]
            if sg(i, w) == 0:
                qv_ret[i] -= 1
            return tuple(permtrim(ret_u)), v, tuple(permtrim(ret_w)), qv_ret, True
    return u, v, w, qv, False


def longest_element(indices):
    perm = [1, 2]
    did_one = True
    while did_one:
        did_one = False
        for i in range(len(indices)):
            j = indices[i] - 1
            if sg(j, perm) == 0:
                if len(perm) < j + 2:
                    perm = perm + list(range(len(perm) + 1, j + 3))
                perm[j], perm[j + 1] = perm[j + 1], perm[j]
                did_one = True
    return permtrim(perm)


def count_less_than(arr, val):
    ct = 0
    i = 0
    while i < len(arr) and arr[i] < val:
        i += 1
        ct += 1
    return ct


def is_parabolic(w, parabolic_index):
    for i in parabolic_index:
        if sg(i - 1, w) == 1:
            return False
    return True


def check_blocks(qv, parabolic_index):
    blocks = []
    cur_block = []
    last_val = -1
    for i in range(len(parabolic_index)):
        if last_val == -1 or last_val + 1 == parabolic_index[i]:
            last_val = parabolic_index[i]
            cur_block += [last_val]
        else:
            blocks += [cur_block]
            cur_block = []
    for block in blocks:
        for i in range(len(block)):
            for j in range(i, len(block)):
                val = 0
                for k in range(i, j + 1):
                    val += omega(block[k], qv)
                if val != 0 and val != -1:
                    return False
    return True


# perms and inversion diff
def kdown_perms(perm, monoperm, p, k):
    inv_m = inv(monoperm)
    inv_p = inv(perm)
    full_perm_list = []
    print(f"LEGACY {perm=} {monoperm=} {inv_m=} {inv_p=} {mulperm(list(perm),monoperm)}")
    if inv(mulperm(list(perm), monoperm)) == inv_m - inv_p:
        full_perm_list += [(tuple(perm), 0, 1)]

    down_perm_list = [(perm, 1)]
    if len(perm) < k:
        return full_perm_list
    a2 = k - 1
    for pp in range(1, p + 1):
        down_perm_list2 = []
        for perm2, s in down_perm_list:
            L = len(perm2)
            print(f"LEGACY {perm2=} {L=}")
            if k > L:
                continue
            s2 = -s
            for b in chain(range(k - 1), range(k, L)):
                if perm2[b] != perm[b]:
                    continue
                if b < a2:
                    i, j = b, a2
                else:
                    i, j, s2 = a2, b, s
                print(f"LEGACY {perm2=} {i=} {j=}")
                if has_bruhat_descent(perm2, i, j):
                    print(f"LEGACY YEAH BABY {perm2=} {i=} {j=}")
                    new_perm = [*perm2]
                    new_perm[a2], new_perm[b] = new_perm[b], new_perm[a2]
                    permtrim(new_perm)
                    print(f"LEGACY {new_perm=}")
                    down_perm_list2 += [(new_perm, s2)]
                    if inv(mulperm(new_perm, monoperm)) == inv_m - inv_p + pp:
                        full_perm_list += [(tuple(new_perm), pp, s2)]
                else:
                    print(f"LEGACY NO BABY {perm2=} {i=} {j=}")
        down_perm_list = down_perm_list2
    return full_perm_list


def compute_vpathdicts(th, vmu, smpify=False):
    vpathdicts = [{} for index in range(len(th))]
    vpathdicts[-1][tuple(vmu)] = None
    thL = len(th)

    top = code(inverse(uncode(th)))
    for i in range(thL - 1, -1, -1):
        top2 = code(inverse(uncode(top)))
        while top2[-1] == 0:
            top2.pop()
        top2.pop()
        top = code(inverse(uncode(top2)))
        monoperm = uncode(top)
        if len(monoperm) < 2:
            monoperm = [1, 2]
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
                if smpify:
                    v2 = sympify(v2)
                vpathdicts2[i][key2].add((key, value[1], v2))
    # print(vpathdicts2)
    return vpathdicts2


def theta(perm):
    cd = code(perm)
    for i in range(len(cd) - 1, 0, -1):
        for j in range(i - 1, -1, -1):
            if cd[j] < cd[i]:
                cd[i] += 1
    cd.sort(reverse=True)
    return cd


def add_perm_dict(d1, d2):
    for k, v in d2.items():
        d1[k] = d1.get(k, 0) + v
    return d1


one = sympify(1)


def elem_sym_poly_q(p, k, varl1, varl2, q_var=q_var):
    if p == 0 and k >= 0:
        return one
    if p < 0 or p > k:
        return zero
    return (
        (varl1[k - 1] - varl2[k - p]) * elem_sym_poly_q(p - 1, k - 1, varl1, varl2, q_var)
        + elem_sym_poly_q(p, k - 1, varl1, varl2, q_var)
        + q_var[k - 1] * elem_sym_poly_q(p - 2, k - 2, varl1, varl2, q_var)
    )


def elem_sym_poly(p, k, varl1, varl2, xstart=0, ystart=0):
    if p > k:
        return zero
    if p == 0:
        return one
    if p == 1:
        res = varl1[xstart] - varl2[ystart]
        for i in range(1, k):
            res += varl1[xstart + i] - varl2[ystart + i]
        return res
    if p == k:
        res = (varl1[xstart] - varl2[ystart]) * (varl1[xstart + 1] - varl2[ystart])
        for i in range(2, k):
            res *= varl1[i + xstart] - varl2[ystart]
        return res
    mid = k // 2
    xsm = xstart + mid
    ysm = ystart + mid
    kmm = k - mid
    res = elem_sym_poly(p, mid, varl1, varl2, xstart, ystart) + elem_sym_poly(
        p,
        kmm,
        varl1,
        varl2,
        xsm,
        ysm,
    )
    for p2 in range(max(1, p - kmm), min(p, mid + 1)):
        res += elem_sym_poly(p2, mid, varl1, varl2, xstart, ystart) * elem_sym_poly(
            p - p2,
            kmm,
            varl1,
            varl2,
            xsm,
            ysm - p2,
        )
    return res


@cache
def call_zvars(v1, v2, k, i):  # noqa: ARG001
    v3 = [*v2, *list(range(len(v2) + 1, i + 1))]
    return [v3[i - 1]] + [v3[j] for j in range(len(v1), len(v3)) if v3[j] != j + 1 and j != i - 1] + [v3[j] for j in range(len(v1)) if v1[j] != v3[j] and j != i - 1]


def elem_sym_func(k, i, u1, u2, v1, v2, udiff, vdiff, varl1, varl2):
    newk = k - udiff
    if newk < vdiff:
        return zero
    if newk == vdiff:
        return one
    yvars = []
    for j in range(min(len(u1), k)):
        if u1[j] == u2[j]:
            yvars += [varl1[u2[j]]]
    for j in range(len(u1), min(k, len(u2))):
        if u2[j] == j + 1:
            yvars += [varl1[u2[j]]]
    for j in range(len(u2), k):
        yvars += [varl1[j + 1]]
    zvars = [varl2[i] for i in call_zvars(v1, v2, k, i)]
    return elem_sym_poly(newk - vdiff, newk, yvars, zvars)


def elem_sym_func_q(k, i, u1, u2, v1, v2, udiff, vdiff, varl1, varl2):
    newk = k - udiff
    if newk < vdiff:
        return zero
    if newk == vdiff:
        return one
    yvars = []
    mlen = max(len(u1), len(u2))
    u1 = [*u1] + [a + 1 for a in range(len(u1), mlen)]
    u2 = [*u2] + [a + 1 for a in range(len(u2), mlen)]
    for j in range(min(len(u1), k)):
        if u1[j] == u2[j]:
            yvars += [varl1[u2[j]]]
    for j in range(len(u1), min(k, len(u2))):
        if u2[j] == j + 1:
            yvars += [varl1[u2[j]]]
    for j in range(len(u2), k):
        yvars += [varl1[j + 1]]
    zvars = [varl2[a] for a in call_zvars(v1, v2, k, i)]
    return elem_sym_poly(newk - vdiff, newk, yvars, zvars)


def trimcode(perm):
    cd = code(perm)
    while len(cd) > 0 and cd[-1] == 0:
        cd.pop()
    return cd


def p_trans(part):
    newpart = []
    if len(part) == 0 or part[0] == 0:
        return [0]
    for i in range(1, part[0] + 1):
        cnt = 0
        for j in range(len(part)):
            if part[j] >= i:
                cnt += 1
        if cnt == 0:
            break
        newpart += [cnt]
    return newpart


def cycle(p, q):
    return list(range(1, p)) + [i + 1 for i in range(p, p + q)] + [p]


def phi1(u):
    c_star = code(inverse(u))
    c_star.pop(0)
    return inverse(uncode(c_star))


def one_dominates(u, w):
    c_star_u = code(inverse(u))
    c_star_w = code(inverse(w))

    a = c_star_u[0]
    b = c_star_w[0]

    for i in range(a, b):
        if i >= len(u) - 1:
            return True
        if u[i] > u[i + 1]:
            return False
    return True


def dominates(u, w):
    u2 = [*u]
    w2 = [*w]
    while u2 != [1, 2] and one_dominates(u2, w2):
        u2 = phi1(u2)
        w2 = phi1(w2)
    if u2 == [1, 2]:
        return True
    return False


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

    if inv(mulperm(list(w), mu_uv_inv)) != inv(mu_uv_inv) - inv(w):
        return u, v, w

    umu = mulperm(list(u), mu_u_inv)
    vmu = mulperm(list(v), mu_v_inv)
    wmu = mulperm(list(w), mu_uv_inv)

    t_mu_w = theta(inverse(wmu))

    mu_w = uncode(t_mu_w)

    w_prime = mulperm(wmu, mu_w)

    if permtrim(list(w)) == permtrim(w_prime):
        return (permtrim(list(u)), permtrim(list(v)), permtrim(list(w)))

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

    return (
        permtrim(mulperm(umu, mu_w_A)),
        permtrim(mulperm(vmu, mu_w_B)),
        permtrim(w_prime),
    )


def mu_A(mu, A):
    mu_t = p_trans(mu)
    mu_A_t = []
    for i in range(len(A)):
        if A[i] < len(mu_t):
            mu_A_t += [mu_t[A[i]]]
    return p_trans(mu_A_t)


def reduce_descents(u, v, w):
    u2 = [*u]
    v2 = [*v]
    w2 = [*w]
    found_one = True
    while found_one:
        found_one = False
        if will_formula_work(u2, v2) or will_formula_work(v2, u2) or one_dominates(u2, w2) or is_reducible(v2) or inv(w2) - inv(u2) == 1:
            break
        for i in range(len(w2) - 2, -1, -1):
            if w2[i] > w2[i + 1] and i < len(v2) - 1 and v2[i] > v2[i + 1] and (i >= len(u2) - 1 or u2[i] < u2[i + 1]):
                w2[i], w2[i + 1] = w2[i + 1], w2[i]
                v2[i], v2[i + 1] = v2[i + 1], v2[i]
                found_one = True
            elif w2[i] > w2[i + 1] and i < len(u2) - 1 and u2[i] > u2[i + 1] and (i >= len(v2) - 1 or v2[i] < v2[i + 1]):
                w2[i], w2[i + 1] = w2[i + 1], w2[i]
                u2[i], u2[i + 1] = u2[i + 1], u2[i]
                found_one = True
            if found_one:
                break
    return permtrim(u2), permtrim(v2), permtrim(w2)


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


def try_reduce_v(u, v, w):
    if is_reducible(v):
        return tuple(permtrim([*u])), tuple(permtrim([*v])), tuple(permtrim([*w]))
    u2 = [*u]
    v2 = [*v]
    w2 = [*w]
    cv = code(v2)
    for i in range(len(v2) - 2, -1, -1):
        if cv[i] == 0 and i < len(cv) - 1 and cv[i + 1] != 0:
            if i >= len(u2) - 1 or u2[i] < u2[i + 1]:
                v2[i], v2[i + 1] = v2[i + 1], v2[i]
                if i >= len(w2) - 1:
                    w2 += list(range(len(w2) + 1, i + 3))
                w2[i + 1], w2[i] = w2[i], w2[i + 1]
                if is_reducible(v2):
                    return tuple(permtrim(u2)), tuple(permtrim(v2)), tuple(permtrim(w2))
                return try_reduce_v(u2, v2, w2)
            if i < len(w2) - 1 and w2[i] > w2[i + 1]:
                u2[i], u2[i + 1] = u2[i + 1], u2[i]
                v2[i], v2[i + 1] = v2[i + 1], v2[i]
                return try_reduce_v(u2, v2, w2)
            return tuple(permtrim(u2)), tuple(permtrim(v2)), tuple(permtrim(w2))
    return tuple(permtrim(u2)), tuple(permtrim(v2)), tuple(permtrim(w2))


def try_reduce_u(u, v, w):
    if one_dominates(u, w):
        return u, v, w
    u2 = [*u]
    v2 = [*v]
    w2 = [*w]
    cu = code(u)
    for i in range(len(u2) - 2, -1, -1):
        if cu[i] == 0 and i < len(cu) - 1 and cu[i + 1] != 0:
            if i >= len(v2) - 1 or v2[i] < v2[i + 1]:
                u2[i], u2[i + 1] = u2[i + 1], u2[i]
                if i > len(w2) - 1:
                    w2 += list(range(len(w2) + 1, i + 3))
                w2[i + 1], w2[i] = w2[i], w2[i + 1]
                if one_dominates(u, w):
                    return tuple(permtrim(u2)), tuple(permtrim(v2)), tuple(permtrim(w2))
                return try_reduce_u(u2, v2, w2)
            if i < len(w2) - 1 and w2[i] > w2[i + 1]:
                u2[i], u2[i + 1] = u2[i + 1], u2[i]
                v2[i], v2[i + 1] = v2[i + 1], v2[i]
                return try_reduce_u(u2, v2, w2)
            return tuple(permtrim(u2)), tuple(permtrim(v2)), tuple(permtrim(w2))
    return tuple(permtrim(u2)), tuple(permtrim(v2)), tuple(permtrim(w2))


def divdiffable(v, u):
    inv_v = inv(v)
    inv_u = inv(u)
    perm2 = permtrim(mulperm(v, inverse(u)))
    if inv(perm2) != inv_v - inv_u:
        return []
    return perm2


def will_formula_work(u, v):
    muv = uncode(theta(v))
    vn1muv = mulperm(inverse(v), muv)
    while True:
        found_one = False
        for i in range(len(vn1muv) - 1):
            if vn1muv[i] > vn1muv[i + 1]:
                found_one = True
                if i < len(u) - 1 and u[i] > u[i + 1]:
                    return False
                vn1muv[i], vn1muv[i + 1] = vn1muv[i + 1], vn1muv[i]
                break
        if not found_one:
            return True


def pull_out_var(vnum, v):
    vup = [*v, len(v) + 1]
    if vnum >= len(v):
        return [[[], v]]
    vpm_list = [(vup, 0)]
    ret_list = []
    for p in range(len(v) + 1 - vnum):
        vpm_list2 = []
        for vpm, b in vpm_list:
            if vpm[vnum - 1] == len(v) + 1:
                vpm2 = [*vpm]
                vpm2.pop(vnum - 1)
                vp = permtrim(vpm2)
                ret_list += [
                    [
                        [v[i] for i in range(vnum, len(v)) if ((i > len(vp) and v[i] == i) or (i <= len(vp) and v[i] == vp[i - 1]))],
                        vp,
                    ],
                ]
            for j in range(vnum, len(vup)):
                if vpm[j] <= b:
                    continue
                for i in range(vnum):
                    if has_bruhat_ascent(vpm, i, j):
                        vpm[i], vpm[j] = vpm[j], vpm[i]
                        vpm_list2 += [([*vpm], vpm[i])]
                        vpm[i], vpm[j] = vpm[j], vpm[i]
        vpm_list = vpm_list2
    for vpm, b in vpm_list:
        if vpm[vnum - 1] == len(v) + 1:
            vpm2 = [*vpm]
            vpm2.pop(vnum - 1)
            vp = permtrim(vpm2)
            ret_list += [
                [
                    [v[i] for i in range(vnum, len(v)) if ((i > len(vp) and v[i] == i) or (i <= len(vp) and v[i] == vp[i - 1]))],
                    vp,
                ],
            ]
    return ret_list


def get_cycles(perm):
    cycle_set = []
    done_vals = set()
    for i in range(len(perm)):
        p = i + 1
        if perm[i] == p:
            continue
        if p in done_vals:
            continue
        cycle = []
        m = -1
        max_index = -1
        while p not in done_vals:
            cycle += [p]
            done_vals.add(p)
            if p > m:
                m = p
                max_index = len(cycle) - 1
            p = perm[p - 1]
        cycle = tuple(cycle[max_index + 1 :] + cycle[: max_index + 1])
        cycle_set += [cycle]
    return cycle_set


def double_elem_sym_q(u, p1, p2, k, q_var=q_var):
    ret_list = {}
    perms1 = elem_sym_perms_q(u, p1, k, q_var)
    iu = inverse(u)
    for perm1, udiff1, mul_val1 in perms1:
        perms2 = elem_sym_perms_q(perm1, p2, k, q_var)
        cycles1 = get_cycles(tuple(permtrim(mulperm(iu, [*perm1]))))
        cycles1_dict = {}
        for c in cycles1:
            if c[-1] not in cycles1_dict:
                cycles1_dict[c[-1]] = []
            cycles1_dict[c[-1]] += [set(c)]
        ip1 = inverse(perm1)
        for perm2, udiff2, mul_val2 in perms2:
            cycles2 = get_cycles(tuple(permtrim(mulperm(ip1, [*perm2]))))
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
                # print(f"{(perm1, udiff1, mul_val1)=}")
                if (perm1, udiff1, mul_val1) not in ret_list:
                    ret_list[(perm1, udiff1, mul_val1)] = []
                ret_list[(perm1, udiff1, mul_val1)] += [(perm2, udiff2, mul_val2)]
    return ret_list


def medium_theta(perm):
    cd = code(perm)
    found_one = True
    while found_one:
        found_one = False
        for i in range(len(cd) - 1):
            if cd[i] < cd[i + 1]:
                found_one = True
                cd[i], cd[i + 1] = cd[i + 1] + 1, cd[i]
                break
            if cd[i] == cd[i + 1] and cd[i] != 0 and i > 0 and cd[i - 1] <= cd[i] + 1:
                # if cd[i]==cd[i+1] and i>0 and cd[i-1]<=cd[i]+1:
                cd[i] += 1
                found_one = True
                break
    return cd
