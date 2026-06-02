"""Count subwords of long_word(n) under two filters and compare."""
from itertools import combinations

from schubmult import Permutation, RCGraph, uncode
from schubmult.combinatorics.double_forest import (
    canonical_forest_from_word,
    forest_from_code,
    long_word,
    sylvester_word,
)
from schubmult.combinatorics.indexed_forests import letterpair, omega_insertion


def _pair_label(values):
    counts = {}
    out = []
    for a in values:
        counts[a] = counts.get(a, 0) + 1
        out.append(letterpair(a, counts[a]))
    return tuple(out)


def _infer_n(code):
    F = forest_from_code(list(code))
    sz = sum(code)
    perm = uncode(list(code))
    from schubmult.combinatorics.double_forest import _max_leaf
    max_leaf = 0
    for tt in F:
        max_leaf = max(max_leaf, _max_leaf(tt))
    return max(max_leaf - 1, sz, len(perm.trimcode), 1)


def count_sylvester_subwords(code, n=None):
    """Count #1: subwords of full (barred+unbarred) long word in Syl(F)."""
    code = list(code)
    F = forest_from_code(code)
    sz = sum(code)
    if n is None:
        n = _infer_n(code)
    target = canonical_forest_from_word(sylvester_word(F))
    word = long_word(n)
    cnt = 0
    for idx in combinations(range(len(word)), sz):
        values = [word[i]["value"] for i in idx]
        if canonical_forest_from_word(values) == target:
            cnt += 1
    return cnt


def count_reduced_with_forest(code, n=None):
    """Count #2: subwords of the FULL long word (barred+unbarred) whose value
    sequence is a reduced word for perm = uncode(code) and whose
    omega_invariant[0] (of reversed values) matches the principal-RC target."""
    code = list(code)
    perm = uncode(code)
    sz = sum(code)
    if n is None:
        n = _infer_n(code)

    word = long_word(n)
    target = omega_insertion(
        _pair_label(tuple(reversed(RCGraph.principal_rc(perm, n).perm_word)))
    )[0]
    cnt = 0
    for idx in combinations(range(len(word)), sz):
        values = tuple(word[i]["value"] for i in idx)
        if Permutation.ref_product(*values) != perm:
            continue
        if omega_insertion(_pair_label(tuple(reversed(values))))[0] != target:
            continue
        cnt += 1
    return cnt


if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    comps = [perm.pad_code(n - 1) for perm in perms]
    all_match = True
    for comp in comps:
        c1 = count_sylvester_subwords(comp)
        c2 = count_reduced_with_forest(comp)
        ok = c1 == c2
        if not ok:
            all_match = False
        print(f"comp={comp}  sylv_count={c1}  refl_count={c2}  match={ok}")
    print(f"\nALL MATCH: {all_match}")
