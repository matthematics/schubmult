from schubmult import *
from schubmult.combinatorics.double_forest import (
    long_word, sylvester_word, forest_from_code, canonical_forest_from_word,
    _max_leaf, letter_weight,
)
from schubmult.combinatorics.indexed_forests import letterpair, omega_insertion
from schubmult.rings.polynomial_algebra import *
from schubmult.symbolic import expand


def _codes_eq_mod_trailing(a, b):
    a = list(a); b = list(b)
    while a and a[-1] == 0:
        a.pop()
    while b and b[-1] == 0:
        b.pop()
    return a == b


def build(code):
    code = list(code)
    perm = uncode(code)
    sz = sum(code)
    F = forest_from_code(code)
    max_leaf = max((_max_leaf(tt) for tt in F), default=0)
    n = max(max_leaf - 1, sz, len(perm.trimcode), 1)
    word = long_word(n)
    target_canon = canonical_forest_from_word(sylvester_word(F))
    return code, perm, sz, F, n, word, target_canon


def subwords_dfs(perm, word):
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


def _pair_label(values):
    counts = {}
    out = []
    for a in values:
        counts[a] = counts.get(a, 0) + 1
        out.append(letterpair(a, counts[a]))
    return tuple(out)


def variant(code, x, t, mode):
    """mode:
       'direct'      : no filter, weight = product of picked letters
       'canon'       : filter canonical_forest_from_word(values)==target, weight=picked letters
       'omega'       : filter omega P-forest code == code, weight=picked letters
    """
    code, perm, sz, F, n, word, target_canon = build(code)
    total = 0
    for subword in subwords_dfs(perm, word):
        values = tuple(lt["value"] for lt in subword)
        if mode == "canon":
            if canonical_forest_from_word(list(values)) != target_canon:
                continue
        elif mode == "omega":
            lbs, dec = omega_insertion(_pair_label(tuple(reversed(values))))
            if not _codes_eq_mod_trailing(lbs.forest.code, code):
                continue
        # weight = product of picked letters' own weights
        w = 1
        for lt in subword:
            w = w * letter_weight(lt, x, t)
        total = total + w
    return expand(total)


if __name__ == "__main__":
    import sys
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    x, t = GeneratingSet("x"), GeneratingSet("t")
    DForest = PolynomialAlgebra(DoubleForestPolyBasis(x, t))
    perms = Permutation.all_permutations(N)
    comps = [perm.pad_code(N - 1) for perm in perms]
    modes = ["direct", "canon", "omega"]
    fails = {m: [] for m in modes}
    for comp in comps:
        forval = DForest(*comp).expand()
        for m in modes:
            try:
                val = variant(comp, x, t, m)
            except Exception as e:
                fails[m].append((tuple(comp), f"EXC:{e}"))
                continue
            if expand(val - forval) != 0:
                fails[m].append((tuple(comp), "mismatch"))
    for m in modes:
        print(f"mode={m}: {len(fails[m])} fails out of {len(comps)}")
        for c, why in fails[m][:12]:
            print(f"    {c}: {why}")
