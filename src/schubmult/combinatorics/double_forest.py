"""Double forest polynomial construction via the vine subword model.

This module provides a non-script home for the core computational routine
`double_forest_polynomial`, suitable for use by library code.
"""

from __future__ import annotations

from itertools import combinations

from sympy import expand


def _rightmost_leaf(t):
    while t[0] == "node":
        t = t[2]
    return t[1]


def _max_leaf(t):
    if t[0] == "leaf":
        return t[1]
    return max(_max_leaf(t[1]), _max_leaf(t[2]))


def _mul_elementary(forest, i):
    """Right-multiply `forest` by elementary forest `i_` in the Thompson monoid."""

    def shift(t):
        if t[0] == "leaf":
            return ("leaf", t[1] + 1) if t[1] > i else t
        return ("node", shift(t[1]), shift(t[2]))

    forest = [shift(t) for t in forest]
    found = [False]

    def wrap(t):
        if found[0]:
            return t
        if t[0] == "leaf":
            if t[1] == i:
                found[0] = True
                return ("node", ("leaf", i), ("leaf", i + 1))
            return t
        nl = wrap(t[1])
        if found[0]:
            return ("node", nl, t[2])
        return ("node", t[1], wrap(t[2]))

    forest = [wrap(t) for t in forest]
    if not found[0]:
        existing = max((_max_leaf(t) for t in forest), default=0)
        for k in range(existing + 1, i):
            forest.append(("leaf", k))
        forest.append(("node", ("leaf", i), ("leaf", i + 1)))
    return forest


def forest_from_code(code):
    """Build indexed forest F with c(F)=code via Thompson monoid factorization."""
    F = []
    for i, c in enumerate(code, start=1):
        for _ in range(c):
            F = _mul_elementary(F, i)
    return tuple(F)


def forest_qdes(code):
    """Return the left terminal set (qdes) of the forest with the given code.

    qdes(F) = {i : c_i > 0 and c_{i+1} = 0}, using 1-based indexing.

    These are the descents of the forest, analogous to descent set of a
    permutation (Nadeau-Spink-Tewari, arXiv:2406.01510, §3.2).
    """
    c = list(code)
    result = []
    for idx, val in enumerate(c):
        if val > 0:
            next_val = c[idx + 1] if idx + 1 < len(c) else 0
            if next_val == 0:
                result.append(idx + 1)  # 1-based
    return result


def forest_code_from_trimming_sequence(trimming_seq):
    """Return the forest composition sfc(F) for the unique indexed forest F
    whose set of trimming sequences contains ``trimming_seq``.

    A trimming sequence (i_1, ..., i_k) for F ∈ IndexedForests is defined
    recursively: i_k ∈ qdes(F) and (i_1,...,i_{k-1}) ∈ Trim(F/i_k).
    (Nadeau-Spink-Tewari, arXiv:2406.01510, Definition 3.8.)

    Recovery uses the inverse operation (blossoming):
        sfc(F · i) = (c_1,...,c_{i-1}, c_i+1, 0, c_{i+1}, c_{i+2},...)
    i.e. increment position i and insert a 0 immediately after.
    Starting from the empty forest and blossoming left-to-right through
    (i_1,...,i_k) reconstructs F.

    Parameters
    ----------
    trimming_seq : iterable of int
        Sequence of positive integers (1-based leaf positions).

    Returns
    -------
    tuple of int
        The composition sfc(F); trailing zeros are kept to preserve the
        ambient length implied by the sequence.

    Examples
    --------
    >>> forest_code_from_trimming_sequence([1, 1, 2, 4, 7])
    (2, 1, 0, 1, 0, 0, 1, 0)
    >>> forest_code_from_trimming_sequence([3])
    (0, 0, 1, 0)
    >>> forest_code_from_trimming_sequence([])
    ()
    """
    sfc = []
    for i in trimming_seq:
        idx = i - 1  # convert to 0-based
        # extend if needed
        while len(sfc) <= idx:
            sfc.append(0)
        # blossom: increment at idx, then insert 0 immediately after
        sfc[idx] += 1
        sfc.insert(idx + 1, 0)
    return tuple(sfc)


def sylvester_word(forest):
    """Return one Sylvester word of `forest` via pre-order traversal."""
    word = []

    def pre(t):
        if t[0] == "leaf":
            return
        word.append(_rightmost_leaf(t[1]))
        pre(t[1])
        pre(t[2])

    for t in forest:
        pre(t)
    return word

def sylvester_forest(code, genset, t):
    from schubmult import RCGraph, uncode
    # forest = forest_from_code(code)
    # word = tuple(reversed(sylvester_word(forest)))
    perm = uncode(code)#Permutation.ref_product(*word)
    return sum([rc0.polyvalue(genset, t) for rc0 in RCGraph.all_rc_graphs(perm, len(code)) if rc0.forest_weight == tuple(code)])

def double_sylvester_forest(code, genset, t):
    import itertools

    from schubmult import DSx, RCGraph, SingleSchubertRing, Sx, uncode
    from schubmult.combinatorics.indexed_forests import letterpair, omega_insertion
    x = genset
    St = SingleSchubertRing(t)
    result = 0
    perm = uncode(code)
    n = len(code)

    target = RCGraph.principal_rc(perm, n).omega_invariant[0]

    def _pair_label(word):
        counts = {}
        out = []
        for a in word:
            aa = int(a)
            counts[aa] = counts.get(aa, 0) + 1
            out.append(letterpair(aa, counts[aa]))
        return tuple(out)

    double_schub1 = Sx([]) * DSx(perm, "t")
    double_schub2 = 0
    for perm0, coeff in double_schub1.items():
        double_schub2 += Sx(perm0) @ St.from_expr(coeff)
    def _joint_weight(rc1, rc2):
        # Vine merged-row weight:
        #   rc1 row i (descending letters): crossing has column = letter - i,
        #     weight (x_{i+1} - t_{letter - i}).
        #   rc2 row i (stored reversed-ascending), read in ascending order: its
        #     k-th crossing (1-indexed) sits to the right of rc1[i] in the merged
        #     diagram with column = len(rc1[i]) + k, weight
        #     (x_{i+1} - t_{len(rc1[i]) + k}).
        w = 1
        for i in range(n):
            row1 = rc1[i]
            row2_asc = list(reversed(rc2[i]))
            offset = len(row1)
            for ell in row1:
                w *= (x[i + 1] - t[ell - i])
            for k, _ell in enumerate(row2_asc, start=1):
                w *= (x[i + 1] - t[offset + k])
        return w

    for (u, v), coeff in double_schub2.items():
        for rc1, rc2 in itertools.product(
            RCGraph.all_rc_graphs(u, n),
            RCGraph.all_rc_graphs(v, n),
        ):
            joint = []
            for i in range(n):
                joint.extend([*rc1[i], *tuple(reversed(rc2[i]))])
            # match the convention of RCGraph.omega_invariant: insert reversed word
            joint_rev = list(reversed(joint))
            inv = omega_insertion(_pair_label(joint_rev))[0]
            if inv != target:
                continue
            result += _joint_weight(rc1, rc2) * coeff
    return result

def _bst_insert(t, a):
    if t is None:
        return (a, None, None)
    r, L, R = t
    if a < r:
        return (r, _bst_insert(L, a), R)
    return (r, L, _bst_insert(R, a))


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


def double_forest_polynomial(code, x_gen, t_gen, n=None):
    """Compute P_F(x;t) for indexed forest F with code c(F)=code."""
    # Preserve the ambient composition length exactly as passed; do not trim
    # trailing zeros (matches ForestPolyBasis conventions).
    code = list(code)
    F = forest_from_code(code)
    target = canonical_forest_from_word(sylvester_word(F))
    sz = sum(code)

    max_leaf = 0
    for t in F:
        max_leaf = max(max_leaf, _max_leaf(t))
    inferred_n = max(max_leaf - 1, sz, 1)
    if n is None:
        n = inferred_n
    elif n < inferred_n:
        raise ValueError(f"n={n} too small; need n >= {inferred_n}")

    word = long_word(n)
    total = 0
    for idx in combinations(range(len(word)), sz):
        values = [word[i]["value"] for i in idx]
        if canonical_forest_from_word(values) != target:
            continue
        w = 1
        for i in idx:
            w = w * letter_weight(word[i], x_gen, t_gen)
        total = total + w
    return expand(total)


def reflection_subword_polynomial(code, x_gen, t_gen, n=None, mode="forest"):
    """Vine subword model from arXiv:2504.15234 (Bergeron-Gagnon-Nadeau-Spink-Tewari).

    Iterates subwords pi of the long word `long_word(n)` of size
        sz = sum(code) = ell(uncode(code))
    and sums wt(pi) under one of two filters on the *value sequence*
    (treating barred and unbarred letters by their value only):

      mode='forest'   :  values must be Sylvester-equivalent to sylvester_word(F),
                         i.e. lie in Syl(F).  Equals P_F(x;t)  (Theorem 5.1).
      mode='schubert' :  values must be a reduced word for perm=uncode(code),
                         i.e. lie in Red(perm).  Equals S_perm(x;t) (Theorem 6.1).

    Weights (Sylvester column convention, paper p.19):
        wt(j^(k))      = x_k - t_j
        wt(barred j^(k)) = t_j - t_k
    realised by `letter_weight` with bracket-indexed `x_gen`, `t_gen`.
    """
    from schubmult import Permutation, uncode

    code = list(code)
    perm = uncode(code)
    sz = sum(code)
    assert sz == perm.inv, f"sum(code)={sz} != perm.inv={perm.inv}"

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

    if mode == "forest":
        sylv_target = canonical_forest_from_word(sylvester_word(F))
        def accept(values):
            return canonical_forest_from_word(list(values)) == sylv_target
    elif mode == "schubert":
        def accept(values):
            return Permutation.ref_product(*values) == perm
    else:
        raise ValueError(f"unknown mode {mode!r}")

    total = 0
    for idx in combinations(range(len(word)), sz):
        values = tuple(word[i]["value"] for i in idx)
        if not accept(values):
            continue
        w = 1
        for i in idx:
            w = w * letter_weight(word[i], x_gen, t_gen)
        total = total + w
    return expand(total)


def reflection_forest_polynomial(code, x_gen, t_gen, n=None, weight_rule=None):
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


def debug_compare_models(code, x_gen, t_gen, n=None):
    """Print, for `code`, the subwords that survive in:
      (A) the alphabet vine model (double_forest_polynomial), and
      (B) the reflection-alphabet model with the omega-insertion filter
          using the principal-RC reversed perm-word as target,
    so we can stare at the symmetric difference.
    """
    from schubmult import Permutation, RCGraph, uncode
    from schubmult.combinatorics.indexed_forests import letterpair, omega_insertion

    code = list(code)
    perm = uncode(code)
    perm_len = perm.inv
    F = forest_from_code(code)
    sz = sum(code)
    max_leaf = 0
    for tt in F:
        max_leaf = max(max_leaf, _max_leaf(tt))
    inferred_n = max(max_leaf - 1, sz, len(perm.trimcode), 1)
    if n is None:
        n = inferred_n

    # ---- Model A: alphabet vine, sylvester filter on values ----
    word_A = long_word(n)
    target_A = canonical_forest_from_word(sylvester_word(F))
    survivors_A = []
    for idx in combinations(range(len(word_A)), sz):
        values = [word_A[i]["value"] for i in idx]
        if canonical_forest_from_word(values) != target_A:
            continue
        weight = 1
        annot = []
        for i in idx:
            lt = word_A[i]
            wt = letter_weight(lt, x_gen, t_gen)
            weight = weight * wt
            annot.append((lt["syllable"], lt["value"], "bar" if lt["barred"] else "x", wt))
        survivors_A.append((tuple(values), annot, expand(weight)))

    # ---- Model B: reflection alphabet, omega-insertion filter ----
    long_w = []
    for k in range(1, n + 1):
        for v in range(n, k - 1, -1):
            long_w.append({"value": v, "syllable": k})

    def _pair_label(values):
        counts = {}
        out = []
        for a in values:
            aa = int(a)
            counts[aa] = counts.get(aa, 0) + 1
            out.append(letterpair(aa, counts[aa]))
        return tuple(out)

    sylv_word_B = tuple(reversed(RCGraph.principal_rc(perm, n).perm_word))
    target_B = omega_insertion(_pair_label(sylv_word_B))[0]
    survivors_B = []
    for idx in combinations(range(len(long_w)), perm_len):
        values = tuple(long_w[i]["value"] for i in idx)
        if Permutation.ref_product(*values) != perm:
            continue
        if omega_insertion(_pair_label(tuple(reversed(values))))[0] != target_B:
            continue
        weight = 1
        annot = []
        for i in idx:
            k = long_w[i]["syllable"]
            v = long_w[i]["value"]
            wt = (x_gen[k] - t_gen[v - k + 1])
            weight = weight * wt
            annot.append((k, v, wt))
        survivors_B.append((tuple(values), annot, expand(weight)))

    print(f"=== code={code} perm={perm} ===")
    print(f"target_A (canonical forest of sylvester_word) = {target_A}")
    print(f"sylv_word_B (reversed principal_rc.perm_word) = {sylv_word_B}")
    print(f"target_B = {target_B}")
    print(f"-- Model A survivors ({len(survivors_A)}):")
    for vals, annot, w in survivors_A:
        print(f"  values={vals}  annot={annot}  weight={w}")
    print(f"-- Model B survivors ({len(survivors_B)}):")
    for vals, annot, w in survivors_B:
        print(f"  values={vals}  annot={annot}  weight={w}")
    sumA = sum((w for _, _, w in survivors_A), 0)
    sumB = sum((w for _, _, w in survivors_B), 0)
    print(f"sum(A)={expand(sumA)}")
    print(f"sum(B)={expand(sumB)}")
    print(f"sum(B) - sum(A) = {expand(sumB - sumA)}")
