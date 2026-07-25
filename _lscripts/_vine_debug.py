from schubmult import *
from schubmult.combinatorics.double_forest import (
    long_word, sylvester_word, forest_from_code, canonical_forest_from_word,
    _max_leaf, letter_weight,
)
from schubmult.combinatorics.indexed_forests import letterpair, omega_insertion
from schubmult.symbolic import expand
from itertools import combinations


def fmt_letter(lt):
    b = "bar" if lt["barred"] else "   "
    return f"{b}{lt['value']}^({lt['syllable']})"


def show_double_forest_model(code, x, t):
    """The KNOWN-correct vine model: subwords of full long word, Sylvester filter on values."""
    code = list(code)
    perm = uncode(code)
    sz = sum(code)
    F = forest_from_code(code)
    target = canonical_forest_from_word(sylvester_word(F))
    max_leaf = max((_max_leaf(tt) for tt in F), default=0)
    n = max(max_leaf - 1, sz, len(perm.trimcode), 1)
    word = long_word(n)
    print(f"=== code={tuple(code)} perm={perm} sz={sz} n={n} ===")
    print(f"sylvester_word(F) = {sylvester_word(F)}  canon target = {target}")
    print("long_word:")
    for i, lt in enumerate(word):
        print(f"   [{i:2}] {fmt_letter(lt)}")
    print("-- correct-model surviving subwords (Sylvester filter on values) --")
    total = 0
    for idx in combinations(range(len(word)), sz):
        values = [word[i]["value"] for i in idx]
        if canonical_forest_from_word(values) != target:
            continue
        w = 1
        picks = []
        for i in idx:
            w = w * letter_weight(word[i], x, t)
            picks.append(fmt_letter(word[i]))
        total = total + w
        print(f"   values={values}  letters={picks}  wt={expand(w)}")
    print(f"   TOTAL = {expand(total)}")
    return expand(total)


def show_reduced_subword_omega(code, x, t):
    """Reduced subwords for perm + omega insertion. Show P/Q forests and codes."""
    code = list(code)
    perm = uncode(code)
    sz = sum(code)
    F = forest_from_code(code)
    max_leaf = max((_max_leaf(tt) for tt in F), default=0)
    n = max(max_leaf - 1, sz, len(perm.trimcode), 1)
    word = long_word(n)
    print("-- reduced subwords for perm + omega insertion --")
    print(f"forest F code = {code}")

    def _pair_label(values):
        counts = {}
        out = []
        for a in values:
            counts[a] = counts.get(a, 0) + 1
            out.append(letterpair(a, counts[a]))
        return tuple(out)

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
        values = tuple(lt["value"] for lt in subword)
        letters = [fmt_letter(lt) for lt in subword]
        lbs, dec = omega_insertion(_pair_label(tuple(reversed(values))))
        pcode = lbs.forest.code
        print(f"   values={values} letters={letters}")
        print(f"       P.code={pcode}  P.support={lbs.forest.support}")
        for nd in lbs.forest.inorder_traversal:
            print(f"         P node idx={nd.index} rho={nd.rho} label={nd.label}")
        for nd, lab in dec.inorder_traversal:
            print(f"         Q node idx={nd.index} rho={nd.rho} Qlabel={lab}")


if __name__ == "__main__":
    x, t = GeneratingSet("x"), GeneratingSet("t")
    for comp in [(0, 1), (1, 0), (0, 2), (1, 1)]:
        show_double_forest_model(comp, x, t)
        show_reduced_subword_omega(comp, x, t)
        print()
