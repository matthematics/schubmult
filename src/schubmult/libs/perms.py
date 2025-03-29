from bisect import bisect_left


# from functools import cache
def getpermval(perm, index):
    if index < len(perm):
        return perm[index]
    return index + 1


def inv(perm):
    L = len(perm)
    if L == 0:
        return 0
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
    ret = [perm1[perm2[i] - 1] for i in range(len(perm2))] + perm1[len(perm2) :]
    # print(f"{ret=}")
    return ret


def uncode(cd):
    cd2 = [*cd]
    if cd2 == []:
        return []
    max_required = max([cd2[i] + i for i in range(len(cd2))])
    cd2 += [0 for i in range(len(cd2), max_required)]
    fullperm = [i + 1 for i in range(len(cd2) + 1)]
    perm = []
    for i in range(len(cd2)):
        perm += [fullperm.pop(cd2[i])]
    perm += [fullperm[0]]
    return perm


def inverse(perm):
    if len(perm) == 0:
        return perm
    retperm = [0 for i in range(len(perm))]
    for i in range(len(perm)):
        retperm[perm[i] - 1] = i + 1
    return retperm


def permtrim(perm):
    L = len(perm)
    while L > 0 and perm[-1] == L:
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
    c_star = trimcode(inverse(u))
    if len(c_star) > 0:
        c_star.pop(0)
    return permtrim(inverse(uncode(c_star)))


def one_dominates(u, w):
    c_star_u = trimcode(inverse(u))
    c_star_w = trimcode(inverse(w))

    a = c_star_u[0]
    b = c_star_w[0]

    for i in range(a, b):
        if i >= len(u) - 1:
            continue
        if u[i] > u[i + 1]:
            return False
    return True


def dominates(u, w):
    u2 = [*u]
    w2 = [*w]
    while inv(u2) > 0 and one_dominates(u2, w2):
        u2 = phi1(u2)
        w2 = phi1(w2)
    if inv(u2) == 0:
        return True
    return False


def mu_A(mu, A):
    mu_t = p_trans(mu)
    mu_A_t = []
    for i in range(len(A)):
        if A[i] < len(mu_t):
            mu_A_t += [mu_t[A[i]]]
    return p_trans(mu_A_t)


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
