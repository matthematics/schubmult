"""Low-level helpers on permutations given as lists/tuples, plus composition and reduced-word utilities.

Most functions here are used by the multiplication kernels in `schubmult.mult` and by
`schubmult.combinatorics.permutation.Permutation`; `add_perm_dict` is the standard way to merge
``{key: coeff}`` expansions.
"""

from __future__ import annotations

from bisect import bisect_left
from collections.abc import Hashable, Iterable, Mapping, Sequence
from typing import TYPE_CHECKING, Any, TypeVar

if TYPE_CHECKING:
    from schubmult.combinatorics.permutation import Permutation

K = TypeVar("K", bound=Hashable)


def permtrim_list(perm: list[int]) -> list[int]:
    """Strip trailing fixed points ``perm[L-1] == L`` from a list in place and return it."""
    L = len(perm)
    while L > 0 and perm[-1] == L:
        L = perm.pop() - 1
    return perm


def has_bruhat_descent(perm: Permutation | Sequence[int], i: int, j: int) -> bool:
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


def count_bruhat(perm: Permutation | Sequence[int], i: int, j: int) -> int:
    """Signed length change ``inv(perm * t_{ij}) - inv(perm)`` for the transposition of positions ``i < j``."""
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


def has_bruhat_ascent(perm: Permutation | Sequence[int], i: int, j: int) -> bool:
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


def omega(i: int, qv: Sequence[int]) -> int:
    """``i``-th entry (1-indexed) of the Cartan-matrix image of the q-exponent vector ``qv``:
    ``2 qv[i] - qv[i-1] - qv[i+1]`` with boundary conventions. Used to convert quantum ``q``
    monomials to weights in the parabolic quantum product.
    """
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


def sg(i: int, w: Permutation | Sequence[int]) -> int:
    """1 if ``w`` has a descent at 0-indexed position ``i``, else 0."""
    if i >= len(w) - 1 or w[i] < w[i + 1]:
        return 0
    return 1


def count_less_than(arr: Sequence[int], val: int) -> int:
    """Number of leading entries of the sorted list ``arr`` that are ``< val``."""
    ct = 0
    i = 0
    while i < len(arr) and arr[i] < val:
        i += 1
        ct += 1
    return ct


def artin_sequences(n: int) -> set[tuple[int, ...]]:
    """All tuples ``(a_1, ..., a_n)`` with ``0 <= a_i <= n + 1 - i`` (Lehmer codes of ``S_{n+1}``)."""
    if n == 0:
        return {()}
    old_seqs = artin_sequences(n - 1)

    ret: set[tuple[int, ...]] = set()
    for seq in old_seqs:
        for i in range(n + 1):
            ret.add((i, *seq))
    return ret


def weak_compositions(length: int, max_degree: int) -> set[tuple[int, ...]]:
    """All tuples of the given ``length`` with entries in ``0..max_degree``."""
    if length == 0:
        return {()}
    old_seqs = weak_compositions(length - 1, max_degree)

    ret: set[tuple[int, ...]] = set()
    for seq in old_seqs:
        for i in range(max_degree + 1):
            ret.add((i, *seq))
    return ret


def is_parabolic(w: Permutation | Sequence[int], parabolic_index: Iterable[int]) -> bool:
    """Whether ``w`` has no descent at any of the (1-indexed) positions in ``parabolic_index``."""
    for i in parabolic_index:
        if sg(i - 1, w) == 1:
            return False
    return True


def add_perm_dict(d1: Mapping[K, Any], d2: Mapping[K, Any]) -> dict[K, Any]:
    """Return ``d1 + d2`` as coefficient dicts (keys merged, values added)."""
    d_ret = {**d1}
    for k, v in d2.items():
        d_ret[k] = d_ret.get(k, 0) + v
    return d_ret


def add_perm_dict_with_coeff(d1: Mapping[K, Any], d2: Mapping[K, Any], coeff: Any) -> dict[K, Any]:
    """Return ``d1 + coeff * d2`` as coefficient dicts."""
    d_ret = {**d1}
    for k, v in d2.items():
        d_ret[k] = d_ret.get(k, 0) + v * coeff
    return d_ret


def p_trans(part: Sequence[int]) -> list[int]:
    """Conjugate (transpose) of a partition given as a weakly decreasing list; ``[0]`` for the empty partition."""
    newpart: list[int] = []
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


def mu_A(mu: Sequence[int], A: Sequence[int]) -> list[int]:
    """The partition whose conjugate consists of the columns of ``mu`` indexed by ``A`` (0-indexed)."""
    mu_t = p_trans(mu)
    mu_A_t: list[int] = []
    for i in range(len(A)):
        if A[i] < len(mu_t):
            mu_A_t += [mu_t[A[i]]]
    return p_trans(mu_A_t)


def get_cycles(perm: Permutation) -> list[tuple[int, ...]]:
    """``perm.get_cycles()``."""
    return perm.get_cycles()


def old_code(perm: Sequence[int]) -> list[int]:
    """Lehmer code of a permutation list computed by successive deletion from ``[1..L]``."""
    L = len(perm)
    ret: list[int] = []
    v = list(range(1, L + 1))
    for i in range(L - 1):
        itr = bisect_left(v, perm[i])
        ret += [itr]
        v = v[:itr] + v[itr + 1 :]
    return ret


def cyclic_sort(L: list[int]) -> list[int]:
    """Rotate the list so its maximum is last."""
    m = max(L)
    i = L.index(m)
    return L[i + 1 :] + L[: i + 1]


def cyclic_sort_min(L: list[int]) -> list[int]:
    """Rotate the list so its minimum is first."""
    m = min(L)
    i = L.index(m)
    return L[i:] + L[:i]


def h_vector(q_vector: Sequence[int]) -> tuple[int, ...]:
    """Positions (1-indexed) where the vector strictly increases, up to its first decrease."""
    h: list[int] = []
    val = 0
    for i in range(len(q_vector)):
        val2 = q_vector[i]
        if val2 < val:
            break
        if val2 > val:
            h += [i + 1]
        val = val2
    return tuple(h)


def l_vector(q_vector: Sequence[int]) -> tuple[int, ...]:
    """Find l_j = last position where d equals j (where d decreases from j to j-1)."""
    l: list[int] = []
    val = 0
    for i in range(len(q_vector)):
        val2 = q_vector[i]
        if val2 < val:
            # Record the PREVIOUS position (i) as the last position with value val
            l += [i]
        val = val2
    return tuple(reversed(l))


def tau_d(d: Sequence[int]) -> Permutation:
    """Partial permutation built from `h_vector`/`l_vector` of ``d`` (``tau[l_i - i] = h_i``), completed by
    ``Permutation.from_partial``.
    """
    from schubmult.combinatorics.permutation import Permutation

    lv = l_vector(d)
    hv = h_vector(d)

    tau: list[int | None] = [None] * len(d)
    for i in range(len(lv)):
        if lv[i] - i >= len(d):
            tau += [None] * (lv[i] - i - len(d) + 1)
        tau[lv[i] - i] = hv[i]
    return Permutation.from_partial(tau)


def phi_d(d: Sequence[int]) -> Permutation:
    """Companion of `tau_d` with the shifted placement ``phi[l_i - 1 - i] = h_i``."""
    from schubmult.combinatorics.permutation import Permutation

    hv = h_vector(d)
    lv = l_vector(d)

    phi: list[int | None] = [None] * len(d)
    for i in range(len(hv)):
        if lv[i] - 1 - i >= len(d):
            phi += [None] * (lv[i] - 1 - i - len(d) + 1)
        phi[lv[i] - 1 - i] = hv[i]
    return Permutation.from_partial(phi)


def conjugate_weak_composition(comp: Sequence[int]) -> tuple[int, ...]:
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


def find_reduced_fail(word: Sequence[int], inserted: int) -> int | None:
    """After changing letter ``inserted`` of a word, find the other position carrying the same root
    (the letter whose deletion would make the word reduced again), or ``None``.
    """
    from schubmult.combinatorics.permutation import Permutation

    perm = Permutation.ref_product(*word)
    a_start, b_start = perm.right_root_at(inserted, word=word)
    # positive = False
    # if a_start > b_start:
    #     positive = True
    # for i in range(len(word)):
    #     if i == inserted:
    #         continue
    #     a, b = perm.right_root_at(i, word=word)
    #     if not positive and a > b:
    #         return i
    #     if positive and a == b_start and b == a_start:
    #         return i
    # return None
    return next(iter([i for i in range(len(word)) if set(perm.right_root_at(i, word=word)) == {a_start, b_start} and i != inserted]), None)


def is_reduced(word: Sequence[int]) -> bool:
    """Whether the word of simple reflections is reduced (``inv`` of its product equals its length)."""
    from schubmult.combinatorics.permutation import Permutation

    return Permutation.ref_product(*word).inv == len(word)


def little_bump_pos(word: Sequence[int], index: int) -> tuple[int, ...]:
    """Little bump at position ``index``: decrement that letter (increment if it is 1), and while the
    word is not reduced, repeat at the letter found by `find_reduced_fail`.
    """

    if not is_reduced(word):
        raise ValueError(f"Word {word} is not reduced, cannot perform Little bump.")
    if index < 0 or index >= len(word):
        raise ValueError(f"Index {index} is out of bounds for word of length {len(word)}.")
    letters = [*word]
    while True:
        if letters[index] == 1:
            letters[index] = letters[index] + 1
        else:
            letters[index] = letters[index] - 1
        if is_reduced(letters):
            break
        fail = find_reduced_fail(letters, index)
        if fail is None:
            raise ValueError(f"Word {letters} is not reduced but no repeated root was found.")
        index = fail
    return tuple(letters)


def little_bump(word: Sequence[int], i: int, j: int) -> tuple[int, ...]:
    """Little bump of a reduced word at the letter whose right root is the inversion ``(i, j)``."""
    from schubmult.combinatorics.permutation import Permutation

    if not is_reduced(word):
        raise ValueError(f"Word {word} is not reduced, cannot perform Little bump.")
    letters = [*word]
    roots = {Permutation._right_root_at(index, letters): index for index in range(len(letters))}
    index = roots.get((i, j), None)
    if index is None:
        raise ValueError(f"Word {letters} does not have an inversion at ({i}, {j})")
    return little_bump_pos(letters, index)


def little_zero(word: Sequence[int], length: int) -> tuple[int, ...]:
    """Repeatedly Little-bump at the last descent until the product's code has fewer than ``length``
    entries (Little's map toward a smaller permutation).
    """
    from schubmult import Permutation

    perm = Permutation.ref_product(*word)
    if len(perm.trimcode) < length:
        return tuple(word)
    if len(perm.trimcode) > length:
        raise ValueError("Word is too long for the specified length")
    new_word: Sequence[int] = [*word]
    while len(perm.trimcode) >= length:
        d = len(perm.trimcode)
        old_word = new_word
        new_word = little_bump(new_word, d, d + 1)
        if new_word == old_word:
            raise ValueError(f"Word cannot be bumped further {word=} {new_word} {length=}")
        perm = Permutation.ref_product(*new_word)
    return tuple(new_word)
