from schubmult import *
from schubmult.combinatorics.indexed_forests import weak_composition_to_indfor, comp_to_thompson_word
from schubmult.rings.polynomial_algebra import *


#(comp, x, y):
def canonical_forest_from_word(word):
    """Canonical forest form used to test Sylvester-equivalence of words."""
    if not word:
        return ()
    sup = sorted(set(word))
    intervals = []
    cur = [sup[0]]
    for x in sup[1:]:
        if x == cur[-1] + 1:
            cur.append(x)
        else:
            intervals.append((cur[0], cur[-1]))
            cur = [x]
    intervals.append((cur[0], cur[-1]))
    out = []
    for lo, hi in intervals:
        iset = set(range(lo, hi + 1))
        sub = [a for a in word if a in iset]
        bst = None
        for a in sub:
            bst = _bst_insert(bst, a)
        out.append((lo, hi, bst))
    return tuple(out)


def long_word(n):
    """Build omega_tilde_[n] as a list of letter records."""
    letters = []
    for k in range(1, n + 1):
        for v in range(n, k - 1, -1):
            letters.append({"value": v, "barred": False, "syllable": k})
        for v in range(k + 1, n + 1):
            letters.append({"value": v, "barred": True, "syllable": k})
    return letters


def letter_weight(letter, x_gen, t_gen):
    """Weight map in the vine model. x_gen, t_gen are 1-indexed subscriptables."""
    i = letter["syllable"]
    j = letter["value"]
    if letter["barred"]:
        return t_gen[j] - t_gen[i]
    return x_gen[i] - t_gen[j]

    
def fake_dforest(code, x_gen, t_gen, n=None, weight_rule=None):
    """Reflection-style forest model from the FULL long word (barred + unbarred).

    Iterates subwords of `long_word(n)` (same enumeration used by
    `double_forest_polynomial`). A subword is kept iff:
      (a) its value sequence is a reduced expression for the same
          permutation as `sylvester_word(F)`, AND
      (b) its omega-insertion P-symbol (Nadeau-Tewari arXiv:2306.10939 §5.1)
          equals that of `sylvester_word(F)` (single Omega-class).

    Default weight is `letter_weight` (paper convention). Override with
    `weight_rule(letter, position_in_subword, picked_letters, x_gen, t_gen,
                 sylvester_letters)` -> sympy expression, where
        letter            = picked long_word letter dict
        position_in_subword = 0-indexed position in the picked subword
        picked_letters    = tuple of all picked long_word letter dicts
        sylvester_letters = tuple of long_word letters of the canonical
                            Sylvester subword that gave `target`
    """
    from schubmult import Permutation, RCGraph, uncode
    from schubmult.combinatorics.indexed_forests import letterpair, omega_insertion

    code = list(code)
    perm = uncode(code)
    sz = sum(code)
    assert sz == perm.inv

    F = forest_from_code(code)
    max_leaf = 0
    for tt in F:
        max_leaf = max(max_leaf, _max_leaf(tt))
    inferred_n = max(max_leaf - 1, sz, len(perm.trimcode), 1)
    if n is None:
        n = inferred_n
    elif n < inferred_n:
        raise ValueError(f"n={n} too small; need n >= {inferred_n}")

    word = long_word(n)

    def _pair_label(values):
        counts = {}
        out = []
        for a in values:
            counts[a] = counts.get(a, 0) + 1
            out.append(letterpair(a, counts[a]))
        return tuple(out)

    sw = sylvester_word(F)
    canonical_red = list(RCGraph.principal_rc(perm, n).perm_word)
    target = omega_insertion(_pair_label(tuple(reversed(canonical_red))))[0]

    sylv_target_canon = canonical_forest_from_word(sw)
    sylvester_letters = None
    for idx in combinations(range(len(word)), sz):
        vs = [word[i]["value"] for i in idx]
        if canonical_forest_from_word(vs) == sylv_target_canon:
            sylvester_letters = tuple(word[i] for i in idx)
            break

    if weight_rule is None:
        def weight_rule(letter, pos, picked, x, t, sylv):  # noqa: ARG001
            return letter_weight(letter, x, t)

    total = 0
    for idx in combinations(range(len(word)), sz):
        values = tuple(word[i]["value"] for i in idx)
        if Permutation.ref_product(*values) != perm:
            continue
        if omega_insertion(_pair_label(tuple(reversed(values))))[0] != target:
            continue
        picked = tuple(word[i] for i in idx)
        w = 1
        for pos, lt in enumerate(picked):
            w = w * weight_rule(lt, pos, picked, x_gen, t_gen, sylvester_letters)
        total = total + w
    return expand(total)

if __name__ == "__main__":
    import sys
    from schubmult.abc import x, y

    n = int(sys.argv[1])
    DForest = PolynomialAlgebra(DoubleForestPolyBasis(x,y))
    subs_dict = {x[i]: y[i] for i in range(100)}
    perms = Permutation.all_permutations(n)
    comps = [perm.pad_code(n - 1) for perm in perms]
    ForestX = PolynomialAlgebra(ForestPolyBasis(x))
    ForestY = PolynomialAlgebra(ForestPolyBasis(y))
    for comp in comps:
        if sum(comp) == 0:
            continue
        poo = fake_dforest(comp, x, y)

        # stinkbat = ForestX.from_expr(poo)
        # new_stinkbat = 0
        # for comp2, coeff in stinkbat.items():
        #     new_stinkbat += ForestX(*comp2) @ ForestY.from_expr(coeff)
        # print(f"{comp}: {new_stinkbat}")

        forest = weak_composition_to_indfor(comp)

        desc = forest.trim_descents[0]
        minus_dict = {x[desc]: 0} | {x[i + 1]: x[i] for i in range(desc, n + 5)}
        plus_dict = {x[desc + 1]: 0} | {x[i + 1]: x[i] for i in range(desc + 1, n + 5)}
        # #minus_dict = {}
        
        minus_dict_y = {y[i + 1]: y[i] for i in range(desc, n + 5)}
        poo_spinach = (((poo.subs(plus_dict).expand() - poo.subs(minus_dict).expand())/(x[desc]))).expand()#.simplify()#.subs(minus_dict_y).expand()
        tforest = forest.trim_descent(desc)
        tcomp = tforest.code
        poo_spinach -= fake_dforest(tcomp, x, y)
        poo_spinach = poo_spinach.expand()
        # #poly = poo.subs({k: v for k, v in subs_dict.items() if k in poo.free_symbols}).expand()
        # #assert poly.as_coefficients_dict() == {comp: 1}
        assert poo_spinach.expand() == 0, f"Error: Double forest polynomial for {comp} is not zero after substitution: {poo_spinach=}\n{poo=}"