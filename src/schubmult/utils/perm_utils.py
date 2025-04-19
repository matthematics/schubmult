from bisect import bisect_left

from schubmult import Permutation


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


@ensure_perms
def code(perm):
    return perm.code


@ensure_perms
def mulperm(perm1, perm2):
    return perm1 * perm2


def uncode(cd):
    cd2 = [*cd]
    if cd2 == []:
        return Permutation([])
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
    while L > 0 and perm[-1] == L:
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


def longest_element(indices):
    perm = Permutation([1, 2])
    did_one = True
    while did_one:
        did_one = False
        for i in range(len(indices)):
            j = indices[i] - 1
            if sg(j, perm) == 0:
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


def mu_A(mu, A):
    mu_t = p_trans(mu)
    mu_A_t = []
    for i in range(len(A)):
        if A[i] < len(mu_t):
            mu_A_t += [mu_t[A[i]]]
    return p_trans(mu_A_t)


def get_cycles(perm):
    return perm.get_cycles()


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


def split_perms(perms):
    perms2 = [perms[0]]
    for perm in perms[1:]:
        cd = code(perm)
        index = -1
        not_zero = False
        did = False
        for i in range(len(cd)):
            if cd[i] != 0:
                not_zero = True
            elif not_zero and cd[i] == 0:
                not_zero = False
                index = i
                num_zeros_to_miss = 0
                for j in range(index):
                    if cd[j] != 0:
                        num_zeros_to_miss = max(num_zeros_to_miss, cd[j] - (index - 1 - j))
                num_zeros = 0
                for j in range(index, len(cd)):
                    if cd[j] != 0:
                        break
                    num_zeros += 1
                if num_zeros >= num_zeros_to_miss:
                    cd1 = cd[:index]
                    cd2 = [0 for i in range(index)] + cd[index:]
                    perms2 += [
                        uncode(cd1),
                        uncode(cd2),
                    ]
                    did = True
                    break
        if not did:
            perms2 += [perm]
    return perms2


def cyclic_sort(L):
    m = max(L)
    i = L.index(m)
    return L[i + 1 :] + L[: i + 1]
