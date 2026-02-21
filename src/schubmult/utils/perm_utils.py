from bisect import bisect_left


def permtrim_list(perm):
    L = len(perm)
    while L > 0 and perm[-1] == L:
        L = perm.pop() - 1
    return perm


def has_bruhat_descent(perm, i, j):
    """Check if perm has a Bruhat descent from position i to j.

    Optimized version assuming perm is a Permutation object with direct indexing.
    """
    perm_i = perm[i]
    perm_j = perm[j]
    if perm_i < perm_j:
        return False
    # Check if there's any value between perm[i] and perm[j] in positions i+1 to j-1
    for p in range(i + 1, j):
        perm_p = perm[p]
        if perm_p > perm_j and perm_i > perm_p:
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
    """Check if perm has a Bruhat ascent from position i to j.

    Optimized version assuming perm is a Permutation object with direct indexing.
    """
    perm_i = perm[i]
    perm_j = perm[j]
    if perm_i > perm_j:
        return False
    # Check if there's any value between perm[i] and perm[j] in positions i+1 to j-1
    for p in range(i + 1, j):
        perm_p = perm[p]
        if perm_i < perm_p < perm_j:
            return False
    return True


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


def count_less_than(arr, val):
    ct = 0
    i = 0
    while i < len(arr) and arr[i] < val:
        i += 1
        ct += 1
    return ct


def artin_sequences(n):
    if n == 0:
        return {()}
    old_seqs = artin_sequences(n - 1)

    ret = set()
    for seq in old_seqs:
        for i in range(n + 1):
            ret.add((i, *seq))
    return ret


def is_parabolic(w, parabolic_index):
    for i in parabolic_index:
        if sg(i - 1, w) == 1:
            return False
    return True


def add_perm_dict(d1, d2):
    d_ret = {**d1}
    for k, v in d2.items():
        d_ret[k] = d_ret.get(k, 0) + v
    return d_ret

def add_perm_dict_with_coeff(d1, d2, coeff):
    d_ret = {**d1}
    for k, v in d2.items():
        d_ret[k] = d_ret.get(k, 0) + v * coeff
    return d_ret

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


def mu_A(mu, A):
    mu_t = p_trans(mu)
    mu_A_t = []
    for i in range(len(A)):
        if A[i] < len(mu_t):
            mu_A_t += [mu_t[A[i]]]
    return p_trans(mu_A_t)


def get_cycles(perm):
    return perm.get_cycles()


def old_code(perm):
    L = len(perm)
    ret = []
    v = list(range(1, L + 1))
    for i in range(L - 1):
        itr = bisect_left(v, perm[i])
        ret += [itr]
        v = v[:itr] + v[itr + 1 :]
    return ret


def cyclic_sort(L):
    m = max(L)
    i = L.index(m)
    return L[i + 1 :] + L[: i + 1]


def cyclic_sort_min(L):
    m = min(L)
    i = L.index(m)
    return L[i:] + L[:i]


def h_vector(q_vector):
    h = []
    val = 0
    for i in range(len(q_vector)):
        val2 = q_vector[i]
        if val2 < val:
            break
        if val2 > val:
            h += [i + 1]
        val = val2
    return tuple(h)


def l_vector(q_vector):
    """Find l_j = last position where d equals j (where d decreases from j to j-1)."""
    l = []
    val = 0
    for i in range(len(q_vector)):
        val2 = q_vector[i]
        if val2 < val:
            # Record the PREVIOUS position (i) as the last position with value val
            l += [i]
        val = val2
    return tuple(reversed(l))


def tau_d(d):
    from schubmult.combinatorics.permutation import Permutation

    lv = l_vector(d)
    hv = h_vector(d)

    tau = [None] * len(d)
    for i in range(len(lv)):
        if lv[i] - i >= len(d):
            tau += [None] * (lv[i] - i - len(d) + 1)
        tau[lv[i] - i] = hv[i]
    return Permutation.from_partial(tau)


def phi_d(d):
    from schubmult.combinatorics.permutation import Permutation

    hv = h_vector(d)
    lv = l_vector(d)

    phi = [None] * len(d)
    for i in range(len(hv)):
        if lv[i] - 1 - i >= len(d):
            phi += [None] * (lv[i] - 1 - i - len(d) + 1)
        phi[lv[i] - 1 - i] = hv[i]
    return Permutation.from_partial(phi)


def conjugate_weak_composition(comp):
    """Compute the conjugate of a weak composition.

    The conjugate of a weak composition α = (α₁, α₂, ..., αₙ) is the weak composition
    β where βⱼ = |{i : αᵢ ≥ j}|, i.e., βⱼ counts how many parts of α are at least j.

    This is equivalent to transposing the Ferrers diagram of the composition.

    Args:
        comp: A sequence (list, tuple) of non-negative integers representing a weak composition.

    Returns:
        A tuple representing the conjugate weak composition.

    Examples:
        >>> conjugate_weak_composition([3, 1, 0, 2])
        (3, 2, 1)
        >>> conjugate_weak_composition([4, 2, 1])
        (3, 2, 1, 1)
        >>> conjugate_weak_composition([])
        ()
        >>> conjugate_weak_composition([0, 0, 0])
        ()
    """
    if not comp:
        return ()

    max_part = max(comp)
    if max_part == 0:
        return ()

    # Count how many parts are >= j for each j from 1 to max_part
    conjugate = []
    for j in range(1, max_part + 1):
        count = sum(1 for part in comp if part >= j)
        conjugate.append(count)

    return tuple(conjugate)


def find_reduced_fail(word, inserted):
    from schubmult.combinatorics.permutation import Permutation

    perm = Permutation.ref_product(*word)
    a_start, b_start = perm.right_root_at(inserted, word=word)
    positive = False
    if a_start > b_start:
        positive = True
    for i in range(len(word)):
        if i == inserted:
            continue
        a, b = perm.right_root_at(i, word=word)
        if not positive and a > b:
            return i
        if positive and a == b_start and b == a_start:
            return i
    return None


def is_reduced(word):
    from schubmult.combinatorics.permutation import Permutation

    return Permutation.ref_product(*word).inv == len(word)


def little_bump(word, i, j):
    """Perform a Little bump on a reduced word at the inversion (i, j)."""
    from schubmult.combinatorics.permutation import Permutation

    if not is_reduced(word):
        raise ValueError(f"Word {word} is not reduced, cannot perform Little bump.")
    word = [*word]
    roots = {Permutation._right_root_at(index, word): index for index in range(len(word))}
    index = roots.get((i, j), None)
    if index is None:
        raise ValueError(f"Word {word} does not have an inversion at ({i}, {j})")
    while True:
        if word[index] == 1:
            word[index] = word[index] + 1
        else:
            word[index] = word[index] - 1
        if is_reduced(word):
            break
        index = find_reduced_fail(word, index)
    return tuple(word)

def little_zero(word, length):
    from schubmult import Permutation

    perm = Permutation.ref_product(*word)
    if len(perm.trimcode) < length:
        return tuple(word)
    if len(perm.trimcode) > length:
        raise ValueError("Word is too long for the specified length")
    new_word = [*word]
    while len(perm.trimcode) >= length:
        d = len(perm.trimcode)
        old_word = new_word
        new_word = little_bump(new_word, d, d + 1)
        if new_word == old_word:
            raise ValueError(f"Word cannot be bumped further {word=} {new_word} {length=}")
        perm = Permutation.ref_product(*new_word)
    return tuple(new_word)
