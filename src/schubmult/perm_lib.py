from bisect import bisect_left
from functools import cache, cached_property
from itertools import chain

import numpy as np
import sympy.combinatorics.permutations as spp
from symengine import Mul, Pow
from sympy import Basic, IndexedBase, sympify

from schubmult.logging import get_logger

logger = get_logger(__name__)

zero = sympify(0)
n = 100

# q_var = symarray("q", n)
q_var = IndexedBase("q")

# TODO: permutations act


def ensure_perms(func):
    def wrapper(*args):
        return func(*[Permutation(arg) if (isinstance(arg, list) or isinstance(arg, tuple)) else arg for arg in args])

    return wrapper


def getpermval(perm, index):
    if index < len(perm):
        return perm[index]
    return index + 1


@ensure_perms
def inv(perm):
    return perm.inv
    # L = len(perm)
    # v = list(range(1, L + 1))
    # ans = 0
    # for i in range(L):
    #     itr = bisect_left(v, perm[i])
    #     ans += itr
    #     v = v[:itr] + v[itr + 1 :]
    # return ans


def code(perm):
    # L = len(perm)
    # ret = []
    # v = list(range(1, L + 1))
    # for i in range(L - 1):
    #     itr = bisect_left(v, perm[i])
    #     ret += [itr]
    #     v = v[:itr] + v[itr + 1 :]
    # return ret
    return perm.code

@ensure_perms
def mulperm(perm1, perm2):
    return perm1 * perm2
    # raise Exception("Don't do this")
    # if len(perm1) < len(perm2):
    #     return [perm1[perm2[i] - 1] if perm2[i] <= len(perm1) else perm2[i] for i in range(len(perm2))]
    # return [perm1[perm2[i] - 1] for i in range(len(perm2))] + perm1[len(perm2) :]


def uncode(cd):
    cd2 = [*cd]
    if cd2 == []:
        return Permutation([1, 2])
    max_required = max([cd2[i] + i for i in range(len(cd2))])
    cd2 += [0 for i in range(len(cd2), max_required)]
    fullperm = [i + 1 for i in range(len(cd2) + 1)]
    perm = []
    for i in range(len(cd2)):
        perm += [fullperm.pop(cd2[i])]
    perm += [fullperm[0]]
    return Permutation(perm)


@ensure_perms
def inverse(perm):
    return ~perm


def permtrim_list(perm):
    L = len(perm)
    while L > 2 and perm[-1] == L:
        L = perm.pop() - 1
    return perm


def permtrim(perm):
    return Permutation(perm)


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


@ensure_perms
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
            pos_list = [i for i in range(k) if up_perm[i] == orig_perm[i]]
            for j in range(min(max(k + 1, len(up_perm) - 1), last_j), k - 1, -1):
                for i in pos_list:
                    ct = count_bruhat(up_perm, i, j)
                    if ct == 1 or ct == 2 * (i - j) + 1:
                        new_perm_add = up_perm.swap(i, j)
                        new_val = val
                        if ct < 0:
                            new_val *= np.prod([q_var[index] for index in range(i + 1, j + 1)])
                        perm_list += [(new_perm_add, new_val, j)]
                        total_list += [(new_perm_add, pp + 1, new_val)]
        up_perm_list = perm_list
    return total_list


@ensure_perms
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
                            new_val *= np.prod([q_var[index] for index in range(i + 1, j + 1)])
                        perm_list += [(new_perm_add, new_val, j)]
                        total_list += [(new_perm_add, pp + 1, new_val)]
        up_perm_list = perm_list
    return total_list


def q_vector(q_exp, q_var=q_var):
    # qvar_list = q_var.tolist()
    ret = []

    if q_exp == 1:
        return ret
    if q_exp.base == q_var:
        i = q_exp.args[1]
        return [0 for j in range(i - 1)] + [1]
    if isinstance(q_exp, Pow):
        qv = q_exp.args[0]
        expon = int(q_exp.args[1])
        i = qv.args[1]
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


def longest_element(indices):
    perm = [1, 2]
    did_one = True
    while did_one:
        did_one = False
        for i in range(len(indices)):
            j = indices[i] - 1
            if sg(j, perm) == 0:
                # if len(perm) < j + 2:
                #     perm = perm + list(range(len(perm) + 1, j + 3))
                # perm[j], perm[j + 1] = perm[j + 1], perm[j]
                perm = perm.swap(j, j + 1)
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
    perm = Permutation(perm)
    monoperm = Permutation(monoperm)
    inv_m = inv(monoperm)
    inv_p = inv(perm)
    full_perm_list = []
    # perm = Permutation(perm)
    if inv(perm * monoperm) == inv_m - inv_p:
        full_perm_list += [(perm, 0, 1)]

    down_perm_list = [(perm, 1)]
    if len(perm) < k:
        return full_perm_list
    a2 = k - 1
    for pp in range(1, p + 1):
        down_perm_list2 = []
        for perm2, s in down_perm_list:
            L = len(perm2)
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
                if has_bruhat_descent(perm2, i, j):
                    new_perm = perm2.swap(a2, b)
                    down_perm_list2 += [(new_perm, s2)]
                    if inv(new_perm * monoperm) == inv_m - inv_p + pp:
                        full_perm_list += [(new_perm, pp, s2)]
        down_perm_list = down_perm_list2
    return full_perm_list


def compute_vpathdicts(th, vmu, smpify=False):
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
        monoperm = Permutation(uncode(top))
        # if len(monoperm) < 2:
        #     monoperm = [1, 2]
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


@ensure_perms
def theta(perm):
    cd = code(perm)
    for i in range(len(cd) - 1, 0, -1):
        for j in range(i - 1, -1, -1):
            if cd[j] < cd[i]:
                cd[i] += 1
    cd.sort(reverse=True)
    return cd


def add_perm_dict(d1, d2):
    d_ret = {**d1}
    for k, v in d2.items():
        d_ret[k] = d_ret.get(k, 0) + v
    return d_ret


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


@ensure_perms
def trimcode(perm):
    cd = perm.code
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
    return Permutation(list(range(1, p)) + [i + 1 for i in range(p, p + q)] + [p])
    # return Permutation.cycle(p, q)


@ensure_perms
def phi1(u):
    c_star = (~u).code
    c_star.pop(0)
    # print(f"{uncode(c_star)=}")
    return ~(uncode(c_star))


@ensure_perms
def one_dominates(u, w):
    c_star_u = (~u).code
    c_star_w = (~w).code

    a = c_star_u[0]
    b = c_star_w[0]

    for i in range(a, b):
        if i >= len(u) - 1:
            return True
        if u[i] > u[i + 1]:
            return False
    return True


def dominates(u, w):
    u2 = u
    w2 = w
    while inv(u2) > 0 and one_dominates(u2, w2):
        u2 = phi1(u2)
        w2 = phi1(w2)
    if inv(u2) == 0:
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


def mu_A(mu, A):
    mu_t = p_trans(mu)
    mu_A_t = []
    for i in range(len(A)):
        if A[i] < len(mu_t):
            mu_A_t += [mu_t[A[i]]]
    return p_trans(mu_A_t)


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


# @ensure_perms
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


# @ensure_perms
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


@ensure_perms
def divdiffable(v, u):
    inv_v = inv(v)
    inv_u = inv(u)
    perm2 = v * (~u)
    if inv(perm2) != inv_v - inv_u:
        return []
    return perm2


# @ensure_perms
def will_formula_work(u, v):
    u, v = Permutation(u), Permutation(v)
    muv = uncode(theta(v))
    vn1muv = (~v) * muv
    while True:
        found_one = False
        for i in range(len(vn1muv) - 1):
            if vn1muv[i] > vn1muv[i + 1]:
                found_one = True
                if i < len(u) - 1 and u[i] > u[i + 1]:
                    return False
                # vn1muv[i], vn1muv[i + 1] = vn1muv[i + 1], vn1muv[i]
                vn1muv = vn1muv.swap(i, i + 1)
                break
        if not found_one:
            return True
    return False


def pull_out_var(vnum, v):
    import sys

    logger.debug(f"pull_out_var {vnum=} {v=} {code(v)=}")
    vup = v
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
            for j in range(vnum, len(vup) + 2):
                if vpm[j] <= b:
                    continue
                for i in range(vnum):
                    if has_bruhat_ascent(vpm, i, j):
                        vpm_list2 += [(vpm.swap(i, j), vpm[j])]
        vpm_list = vpm_list2
    for vpm, b in vpm_list:
        logger.debug(f"{vpm=} {b=}")
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
    logger.debug(f"{ret_list=}")
    return ret_list


@ensure_perms
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
                if (perm1, udiff1, mul_val1) not in ret_list:
                    ret_list[(perm1, udiff1, mul_val1)] = []
                ret_list[(perm1, udiff1, mul_val1)] += [(perm2, udiff2, mul_val2)]
    return ret_list


@ensure_perms
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


def old_code(perm):
    L = len(perm)
    ret = []
    v = list(range(1, L + 1))
    for i in range(L - 1):
        itr = bisect_left(v, perm[i])
        ret += [itr]
        v = v[:itr] + v[itr + 1 :]
    return ret


# from typing import NamedTuple


class Permutation(Basic):
    def __new__(cls, perm, action=None, sperm=None):
        return Permutation.__xnew_cached__(cls, tuple(perm), action, sperm)

    @staticmethod
    @cache
    def __xnew_cached__(_class, perm, action, sperm):
        return Permutation.__xnew__(_class, perm, action, sperm)

    @staticmethod
    def __xnew__(_class, perm, action, sperm):
        obj = Basic.__new__(_class)
        if isinstance(perm, Permutation):
            # print("this is happening")
            obj._perm = perm._perm
            obj._sperm = perm._sperm
            if action is None:
                obj._action = perm._action
            else:
                obj._action = action
        else:
            p = tuple(permtrim_list([*perm]))
            obj._perm = p
            if len(obj._perm) < 2:
                obj._perm = (1, 2)
            if sperm:
                obj._sperm = sperm
            else:
                obj._sperm = spp.Permutation._af_new([i - 1 for i in p])
            obj._action = action
        return obj

    @property
    def code(self):
        return list(self.cached_code())
    
    @cache
    def cached_code(self):
        return self._sperm.inversion_vector()

    @cached_property
    def inv(self):
        return self._sperm.inversions()

    def swap(self, i, j):
        import sys

        new_perm = [*self._perm]
        # print(f"{new_perm=}",file=sys.stderr)
        if i > j:
            i, j = j, i
        # print(f"OLD {new_perm=} {new_perm[i]=} {new_perm[j]=} {i=} {j=} FWIPO",file=sys.stderr)

        if j >= len(new_perm):
            new_perm += list(range(len(new_perm) + 1, j + 2))
            # print(f"bugs {new_perm=}", file=sys.stderr)
        new_perm[i], new_perm[j] = new_perm[j], new_perm[i]
        # print(f"NEW {new_perm=} {new_perm[i]=} {new_perm[j]=} FWoPO",file=sys.stderr)
        return Permutation(new_perm)

    def __getitem__(self, i):
        if isinstance(i, slice):
            return [self[ii] for ii in range(*i.indices(len(self._perm)))]
        if i >= len(self._perm):
            return i + 1
        return self._perm[i]

    def __setitem__(self, i, v):
        raise NotImplementedError

    def __hash__(self):
        return hash(self._perm)

    def __mul__(self, other):
        # print("yay")
        new_sperm = other._sperm * self._sperm
        new_perm = permtrim_list([new_sperm(i) + 1 for i in range(new_sperm.size)])
        # print(f"{new_perm=}")
        if len(new_perm) != new_sperm.size:
            new_sperm = spp.Permutation._af_new([i - 1 for i in new_perm])
        return Permutation(new_perm, new_sperm)

    def __iter__(self):
        yield from self._perm.__iter__()

    def __getslice__(self, i, j):
        return self._perm[i:j]

    def __str__(self):
        # print("yay")
        return str(self._perm)

    def __add__(self, other):
        # print("yay")
        if not isinstance(other, list):
            raise NotImplementedError
        permlist = [*self._perm, *other]
        try:
            return Permutation(permlist)
        except Exception:
            return permlist

    def __radd__(self, other):
        # print("yay")
        if not isinstance(other, list):
            raise NotImplementedError
        permlist = [*other, *self._perm]
        try:
            return Permutation(permlist)
        except Exception:
            return permlist

    def __eq__(self, other):
        # print("yay")
        if isinstance(other, Permutation):
            return other._perm == self._perm
        if isinstance(other, list):
            return [*self._perm] == other
        if isinstance(other, tuple):
            return self._perm == other
        return False

    def __len__(self):
        # print("yay")
        return len(self._perm)

    def __invert__(self):
        new_sperm = ~(self._sperm)
        new_perm = [new_sperm(i) + 1 for i in range(new_sperm.size)]
        return Permutation(new_perm, new_sperm)

    def __repr__(self):
        return self.__str__()

    # we want to format permuations
    def _sympystr(self, p):
        return f"moo{str(self._perm)}"
