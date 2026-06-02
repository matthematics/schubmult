"""Compare per-subword weights between alphabet-vine (canonical_forest filter)
and reflection model (reduced+omega filter) on the same long_word index sets.

Both models select subwords of long_word(n). For each kept subword, print
(idx, values, per-letter weights, total weight) under each model.
"""
from itertools import combinations

from sympy import expand

from schubmult import GeneratingSet, Permutation, RCGraph, uncode
from schubmult.combinatorics.double_forest import (
    canonical_forest_from_word,
    forest_from_code,
    letter_weight,
    long_word,
    sylvester_word,
    _max_leaf,
)
from schubmult.combinatorics.indexed_forests import letterpair, omega_insertion


def _pair_label(values):
    counts = {}
    out = []
    for a in values:
        counts[a] = counts.get(a, 0) + 1
        out.append(letterpair(a, counts[a]))
    return tuple(out)


def diagnose(code, t_index_unbarred, force_n=None):
    code = list(code)
    perm = uncode(code)
    sz = sum(code)
    F = forest_from_code(code)
    max_leaf = 0
    for tt in F:
        max_leaf = max(max_leaf, _max_leaf(tt))
    n = force_n if force_n is not None else max(max_leaf - 1, sz, len(perm.trimcode), 1)

    x = GeneratingSet("x")
    t = GeneratingSet("t")
    word = long_word(n)

    target_A = canonical_forest_from_word(sylvester_word(F))
    target_B = omega_insertion(
        _pair_label(tuple(reversed(RCGraph.principal_rc(perm, n).perm_word)))
    )[0]

    print(f"=== code={code}  perm={perm}  n={n}  sz={sz} ===")
    print(f"long_word values: {[lt['value'] for lt in word]}")
    print(f"long_word barred: {[lt['barred'] for lt in word]}")
    print(f"long_word syll  : {[lt['syllable'] for lt in word]}")

    keepA = []
    keepB = []
    for idx in combinations(range(len(word)), sz):
        values = tuple(word[i]["value"] for i in idx)
        in_A = canonical_forest_from_word(list(values)) == target_A
        in_B = (
            Permutation.ref_product(*values) == perm
            and omega_insertion(_pair_label(tuple(reversed(values))))[0] == target_B
        )
        if in_A:
            keepA.append(idx)
        if in_B:
            keepB.append(idx)

    setA = set(keepA)
    setB = set(keepB)
    print(f"|A|={len(setA)} |B|={len(setB)} A==B: {setA == setB}")
    if setA != setB:
        print(f"  A\\B = {sorted(setA - setB)}")
        print(f"  B\\A = {sorted(setB - setA)}")

    def weight_A(idx):
        # Alphabet-vine: every letter via letter_weight (paper p.19 wts)
        terms = []
        prod = 1
        for i in idx:
            lt = word[i]
            wt = letter_weight(lt, x, t)
            terms.append((lt["syllable"], lt["value"], lt["barred"], wt))
            prod = prod * wt
        return terms, expand(prod)

    def weight_B(idx):
        terms = []
        prod = 1
        row_count = {}
        for i in idx:
            lt = word[i]
            if not lt["barred"]:
                row_count[lt["syllable"]] = row_count.get(lt["syllable"], 0) + 1
        seen_row = {}
        for i in idx:
            lt = word[i]
            if lt["barred"]:
                wt = letter_weight(lt, x, t)
            else:
                k = lt["syllable"]
                v = lt["value"]
                seen_row[k] = seen_row.get(k, 0) + 1
                j = seen_row[k]
                cnt = row_count[k]
                col = v + cnt - j
                wt = x[k] - t[t_index_unbarred(k, v, col)]
            terms.append((lt["syllable"], lt["value"], lt["barred"], wt))
            prod = prod * wt
        return terms, expand(prod)

    sumA = 0
    sumB = 0
    union = sorted(setA | setB)
    print(f"\n  per-subword comparison (only kept ones):")
    for idx in union:
        in_A = idx in setA
        in_B = idx in setB
        tag = ("A" if in_A else " ") + ("B" if in_B else " ")
        if in_A:
            tA, pA = weight_A(idx)
        else:
            tA, pA = None, 0
        if in_B:
            tB, pB = weight_B(idx)
        else:
            tB, pB = None, 0
        if in_A:
            sumA = sumA + pA
        if in_B:
            sumB = sumB + pB
        vals = tuple(word[i]["value"] for i in idx)
        print(f"  [{tag}] idx={idx} values={vals}")
        if in_A:
            print(f"      A terms: {tA}")
            print(f"      A prod : {pA}")
        if in_B:
            print(f"      B terms: {tB}")
            print(f"      B prod : {pB}")
        if in_A and in_B:
            print(f"      diff B-A: {expand(pB - pA)}")

    print(f"\n  sumA = {expand(sumA)}")
    print(f"  sumB = {expand(sumB)}")
    print(f"  sumB - sumA = {expand(sumB - sumA)}")


if __name__ == "__main__":
    import sys
    code = tuple(int(x) for x in sys.argv[1].split(","))
    force_n = int(sys.argv[2]) if len(sys.argv) > 2 else None
    diagnose(code, t_index_unbarred=lambda k, v, p: p, force_n=force_n)
