from schubmult import *
from schubmult.combinatorics.double_forest import (
    long_word, sylvester_word, forest_from_code, canonical_forest_from_word,
    _max_leaf, letter_weight,
)
from schubmult.combinatorics.indexed_forests import letterpair, omega_insertion
from schubmult.rings.polynomial_algebra import *
from schubmult.symbolic import expand
from itertools import combinations


def fmt(lt):
    b = "b" if lt["barred"] else " "
    return f"{b}{lt['value']}^{lt['syllable']}"


def _pair_label(values):
    counts = {}
    out = []
    for a in values:
        counts[a] = counts.get(a, 0) + 1
        out.append(letterpair(a, counts[a]))
    return tuple(out)


def analyze(code, x, t):
    code = list(code)
    perm = uncode(code)
    sz = sum(code)
    F = forest_from_code(code)
    max_leaf = max((_max_leaf(tt) for tt in F), default=0)
    n = max(max_leaf - 1, sz, len(perm.trimcode), 1)
    word = long_word(n)
    target = canonical_forest_from_word(sylvester_word(F))

    print(f"===== code={tuple(code)} perm={perm} sz={sz} n={n} =====")
    print(f"sylvester_word(F)={sylvester_word(F)}")

    print("--- CORRECT forest-model subwords (Sylvester filter on values) ---")
    correct = []
    for idx in combinations(range(len(word)), sz):
        vals = [word[i]["value"] for i in idx]
        if canonical_forest_from_word(vals) != target:
            continue
        picks = [word[i] for i in idx]
        w = 1
        for p in picks:
            w = w * letter_weight(p, x, t)
        correct.append((tuple(vals), [fmt(p) for p in picks], expand(w)))
    for vals, picks, w in correct:
        print(f"   vals={vals} picks={picks} wt={w}")

    print("--- DFS reduced subwords + omega insertion ---")

    def subwords():
        stack = [(perm, [], word)]
        while stack:
            wp, wsub, rem = stack.pop()
            if wp.inv == 0:
                yield wsub
                continue
            if len(rem) < wp.inv:
                continue
            d = rem[-1]["value"]
            nrem = rem[:-1]
            if wp[d - 1] > wp[d]:
                stack.append((wp.swap(d - 1, d), [rem[-1]] + wsub, nrem))
            stack.append((wp, wsub, nrem))

    for subword in subwords():
        vals = tuple(lt["value"] for lt in subword)
        picks = [fmt(lt) for lt in subword]
        lbs, dec = omega_insertion(_pair_label(tuple(reversed(vals))))
        pcode = lbs.forest.code
        pnodes = [(nd.index, nd.rho, str(nd.label)) for nd in lbs.forest.inorder_traversal]
        qnodes = [(nd.index, nd.rho, lab) for nd, lab in dec.inorder_traversal]
        print(f"   vals={vals} picks={picks}")
        print(f"       Pcode={pcode} match={canonical_forest_from_word(list(vals))==target}")
        print(f"       Pnodes(idx,rho,label)={pnodes}")
        print(f"       Qnodes(idx,rho,Qlabel)={qnodes}")


if __name__ == "__main__":
    x, t = GeneratingSet("x"), GeneratingSet("t")
    for comp in [(2, 1, 0), (0, 2, 1)]:
        analyze(comp, x, t)
        print()
