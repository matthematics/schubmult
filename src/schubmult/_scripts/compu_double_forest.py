from collections import Counter
from itertools import combinations

from schubmult import *
from schubmult.combinatorics.double_forest import (
    _max_leaf,
    canonical_forest_from_word,
    forest_from_code,
    letter_weight,
    long_word,
    reflection_forest_polynomial,
    sylvester_word,
)
from schubmult.combinatorics.indexed_forests import letterpair, omega_insertion
from schubmult.rings.polynomial_algebra import *
from schubmult.rings.schubert.schubert_ring import SingleSchubertRing

x = GeneratingSet("x")
t = GeneratingSet("t")
St = SingleSchubertRing(t)

DForest = PolynomialAlgebra(DoubleForestPolyBasis(x,t))


def _pair_label(values):
    counts = {}
    out = []
    for a in values:
        counts[a] = counts.get(a, 0) + 1
        out.append(letterpair(a, counts[a]))
    return tuple(out)


def row_weight(picked, n):
    cnt = Counter(p["syllable"] for p in picked)
    return tuple(cnt.get(i, 0) for i in range(1, n + 1))


def class_key(values):
    """Omega-class invariant: P-symbol of omega_insertion."""
    return omega_insertion(_pair_label(values))[0]


def build_match(comp, n):
    """Match B-subwords to A-subwords by omega_insertion Q-symbol.

    A := Sylvester subwords of long_word (canonical_forest == sylv_canon).
    B := reduced-word subwords for uncode(comp) in the omega-class of the
         canonical reduced word of perm.

    For each subword W (A or B), compute (P,Q) = omega_insertion(reverse(W)).
    Both A and B subwords have shape F. Two subwords with identical Q place
    their i-th reversed letter at the same node, so a_rev[i] <-> b_rev[i]
    pairs them position-by-position; equivalently a[pos] <-> b[pos] in
    original (long_word) order. The Q-bijection is unique and total.
    """
    F = forest_from_code(list(comp))
    sw = sylvester_word(F)
    perm = uncode(list(comp))
    sylv_canon = canonical_forest_from_word(sw)
    canonical_red = list(RCGraph.principal_rc(perm, n).perm_word)
    target_P = class_key(tuple(reversed(canonical_red)))

    word = long_word(n)
    sz = sum(comp)

    def q_key(values):
        return omega_insertion(_pair_label(tuple(reversed(values))))[1]

    A_by_Q = {}
    B_by_Q = {}
    for idx in combinations(range(len(word)), sz):
        picked = tuple(word[i] for i in idx)
        vs = tuple(p["value"] for p in picked)
        if canonical_forest_from_word(list(vs)) == sylv_canon:
            A_by_Q.setdefault(q_key(vs), []).append((idx, picked))
        if Permutation.ref_product(*vs) == perm and class_key(tuple(reversed(vs))) == target_P:
            B_by_Q.setdefault(q_key(vs), []).append((idx, picked))

    matching = {}
    collisions = []
    for Q, B_list in B_by_Q.items():
        A_list = A_by_Q.get(Q, [])
        if len(A_list) != len(B_list):
            collisions.append(("Q_count_mismatch", len(A_list), len(B_list)))
        for (idxA, pickedA), (idxB, pickedB) in zip(sorted(A_list), sorted(B_list)):
            matching[idxB] = pickedA
    for Q in A_by_Q.keys() - B_by_Q.keys():
        collisions.append(("A_only_Q", len(A_by_Q[Q])))

    return matching, collisions


def make_match_rule(comp, n):
    matching, collisions = build_match(comp, n)
    if collisions:
        print(f"  [collisions for {comp}]: {collisions[:3]}{'...' if len(collisions)>3 else ''}")

    def rule(letter, pos, picked, x_gen, t_gen, sylvester_letters):
        # Identify B's idx-tuple from picked: each picked has unique pos in long_word;
        # but we need the global idx. Reconstruct by looking up picked tuple in matching.
        # Cache lookup via id(picked) wouldn't survive across calls; rebuild key.
        # picked here is the FULL B-subword for every call.
        key = tuple(_long_word_pos(p, comp, n) for p in picked)
        A_picked = matching.get(key)
        if A_picked is None:
            # fallback to literal
            return letter_weight(letter, x_gen, t_gen)
        return letter_weight(A_picked[pos], x_gen, t_gen)
    return rule


_LW_CACHE = {}
def _long_word_pos(letter, comp, n):
    """Find index of `letter` (as object) in long_word(n)."""
    key = n
    if key not in _LW_CACHE:
        _LW_CACHE[key] = long_word(n)
    lw = _LW_CACHE[key]
    for i, l in enumerate(lw):
        if l is letter:
            return i
    # fall back to value/syllable/barred match
    for i, l in enumerate(lw):
        if l == letter:
            return i
    raise KeyError("letter not found")


if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)
    comps = [perm.pad_code(n - 1) for perm in perms]

    def _infer_n(comp):
        F = forest_from_code(list(comp))
        max_leaf = max((_max_leaf(tt) for tt in F), default=0)
        p = uncode(list(comp))
        return max(max_leaf - 1, sum(comp), len(p.trimcode), 1)

    for comp in comps:
        if sum(comp) == 0:
            continue
        to_compute = DForest(*comp).change_basis(MonomialBasis)
        target = to_compute.expand()
        perm = uncode(list(comp))

        nn = _infer_n(comp)
        rule = make_match_rule(comp, nn)
        refl_forest = reflection_forest_polynomial(comp, x, t, n=nn, weight_rule=rule).expand()

        zero_t = {t[i]: 0 for i in range(1, 20)}
        diff_full = (refl_forest - target).expand()
        x_part_refl = refl_forest.subs(zero_t).expand() if hasattr(refl_forest, "subs") else refl_forest
        x_part_target = target.subs(zero_t).expand() if hasattr(target, "subs") else target
        x_ok = (x_part_refl - x_part_target).expand() == 0
        full_ok = diff_full == 0
        print(f"comp={comp} perm={perm} x_match={x_ok} full_match={full_ok}")
        if not full_ok:
            print(f"  refl - forest = {diff_full}")
