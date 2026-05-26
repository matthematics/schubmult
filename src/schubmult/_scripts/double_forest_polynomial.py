"""Compute double forest polynomials P_F(x; t) via the vine subword model.

Reference
---------
N. Bergeron, L. Gagnon, P. Nadeau, H. Spink, V. Tewari,
"Equivariant quasisymmetry and noncrossing partitions" (arXiv:2504.15234),
Theorem 5.1 and Section 5.1.

Given a code c = (c_1, c_2, ...) for an indexed forest F (so that c(F) = c),
and two sequences of symbols x = (x_1, x_2, ...) and t = (t_1, t_2, ...),
this returns the double forest polynomial

    P_F(x; t) = sum_{pi in R(omega_tilde_[n]; F)} wt(pi),

where omega_tilde_[n] is the long vine word and the sum is over its
subwords whose value-sequence is a Sylvester word for F.
"""

from __future__ import annotations

from itertools import combinations

from sympy import IndexedBase, expand, factor


# --------------------------------------------------------------------------- #
# Forest representation                                                       #
# --------------------------------------------------------------------------- #
# A binary tree is either
#   ("leaf", n)        -- a leaf with absolute label n in N, or
#   ("node", L, R)     -- an internal node with binary-tree children L, R.
# A forest is a tuple of binary trees, ordered by leaf positions.


def _rightmost_leaf(t):
    while t[0] == "node":
        t = t[2]
    return t[1]


def _max_leaf(t):
    if t[0] == "leaf":
        return t[1]
    return max(_max_leaf(t[1]), _max_leaf(t[2]))


def _mul_elementary(forest, i):
    """Right-multiply ``forest`` by the elementary forest ``i_`` in the
    Thompson monoid: leaves > i shift up by 1, and leaf i is wrapped by a
    new internal node whose right child is a brand-new leaf i+1."""

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
    """Build the indexed forest F with c(F) = ``code`` via the Thompson
    monoid factorization F = 1^{c_1} . 2^{c_2} . 3^{c_3} . ..."""
    F = []
    for i, c in enumerate(code, start=1):
        for _ in range(c):
            F = _mul_elementary(F, i)
    return tuple(F)


def sylvester_word(forest):
    """Return one Sylvester word of ``forest`` via pre-order traversal
    (root, left subtree, right subtree) over each non-trivial tree."""
    word = []

    def pre(t):
        if t[0] == "leaf":
            return
        word.append(_rightmost_leaf(t[1]))  # canonical label
        pre(t[1])
        pre(t[2])

    for t in forest:
        pre(t)
    return word


# --------------------------------------------------------------------------- #
# Canonical form for testing whether two words are Sylvester-equivalent       #
# --------------------------------------------------------------------------- #
# Every injective word w is a Sylvester word of a unique indexed forest,
# obtained as follows (Section 3.1 of the paper):
#   - Partition supp(w) into maximal contiguous integer intervals.
#   - For each interval I, take the subword of w using only letters in I.
#   - The BST built by inserting that subword left-to-right is the labeled
#     internal-node skeleton of one tree of the forest.
# Two words are Sylvester-equivalent iff their intervals and per-interval
# BSTs coincide.


def _bst_insert(t, a):
    if t is None:
        return (a, None, None)
    r, L, R = t
    if a < r:
        return (r, _bst_insert(L, a), R)
    return (r, L, _bst_insert(R, a))


def canonical_forest_from_word(word):
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


# --------------------------------------------------------------------------- #
# The vine long word and weights                                              #
# --------------------------------------------------------------------------- #


def long_word(n):
    """Build omega_tilde_[n] as a list of letter records.

    omega^(k)_[n] = (n, n-1, ..., k+1, k, k+1_bar, ..., n-1_bar, n_bar)
    Concatenate for k = 1, ..., n.
    """
    letters = []
    for k in range(1, n + 1):
        for v in range(n, k - 1, -1):
            letters.append({"value": v, "barred": False, "syllable": k})
        for v in range(k + 1, n + 1):
            letters.append({"value": v, "barred": True, "syllable": k})
    return letters


def letter_weight(letter, x_gen, t_gen):
    """wt(j^(i))      = x_i - t_j   (unbarred)
    wt(j_bar^(i)) = t_j - t_i   (barred)"""
    i = letter["syllable"]
    j = letter["value"]
    if letter["barred"]:
        return t_gen(j) - t_gen(i)
    return x_gen(i) - t_gen(j)


# --------------------------------------------------------------------------- #
# Main entry point                                                            #
# --------------------------------------------------------------------------- #


def double_forest_polynomial(code, x_gen, t_gen, n=None):
    """Compute P_F(x; t) for the indexed forest F with c(F) = ``code``.

    Parameters
    ----------
    code : sequence[int]
        The code (c_1, c_2, ...) of F. Trailing zeros are ignored.
    x_gen : callable[int -> Expr]
        Maps i -> x_i  (the non-equivariant variables).
    t_gen : callable[int -> Expr]
        Maps j -> t_j  (the equivariant variables).
    n : int, optional
        Use the long word omega_tilde_[n]. If ``None``, the smallest n with
        F supported in {1, ..., n+1} is used.

    Returns
    -------
    sympy.Expr
        The expanded polynomial P_F(x; t).
    """
    code = list(code)
    while code and code[-1] == 0:
        code.pop()
    F = forest_from_code(code)
    target = canonical_forest_from_word(sylvester_word(F))
    sz = sum(code)  # number of internal nodes = |F|

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


# --------------------------------------------------------------------------- #
# Self-tests against examples from the paper                                  #
# --------------------------------------------------------------------------- #


def _selftest():
    x = IndexedBase("x")
    t = IndexedBase("t")
    xg = lambda i: x[i]
    tg = lambda j: t[j]

    # Example 4.4: P_2 = x_1 + x_2 - t_1 - t_2
    p2 = double_forest_polynomial([0, 1], xg, tg)
    expected = x[1] + x[2] - t[1] - t[2]
    assert expand(p2 - expected) == 0, f"P_2 mismatch: {p2}"
    print(f"P_2                  = {p2}")

    # Figure 1: F has code (2, 0, 1)
    # P_F = (x_1 - t_1)(x_1 - t_2)(x_2 + x_3 - t_1 - t_2)
    pF = double_forest_polynomial([2, 0, 1], xg, tg)
    expected = expand((x[1] - t[1]) * (x[1] - t[2]) * (x[2] + x[3] - t[1] - t[2]))
    assert expand(pF - expected) == 0, f"P_(2,0,1) mismatch: {pF}"
    print(f"P_(code=(2,0,1))     = {factor(pF)}")

    # Example 4.11: F has code (0, 2, 1) (a zigzag forest)
    # F_(0,2,1) = (x_1-t_4)(x_1-t_1)(x_2+x_3-t_1-t_2)
    #          + (x_1+x_2-t_1-t_4)(x_2-t_2)(x_3-t_2)
    p = double_forest_polynomial([0, 2, 1], xg, tg)
    expected = expand(
        (x[1] - t[4]) * (x[1] - t[1]) * (x[2] + x[3] - t[1] - t[2])
        + (x[1] + x[2] - t[1] - t[4]) * (x[2] - t[2]) * (x[3] - t[2])
    )
    assert expand(p - expected) == 0, f"P_(0,2,1) mismatch: {p}"
    print(f"P_(code=(0,2,1))     = {p}")

    # Spot-check Table 1: P_3 = (x_1-t_3)+(t_3-t_1)+(x_2-t_3)+(t_3-t_2)+(x_3-t_3)
    p3 = double_forest_polynomial([0, 0, 1], xg, tg)
    expected = (
        (x[1] - t[3]) + (t[3] - t[1]) + (x[2] - t[3]) + (t[3] - t[2]) + (x[3] - t[3])
    )
    assert expand(p3 - expected) == 0
    print(f"P_3                  = {p3}")

    # Empty forest: P_empty = 1
    pe = double_forest_polynomial([], xg, tg)
    assert pe == 1
    print(f"P_empty              = {pe}")

    print("\nAll self-tests passed.")


if __name__ == "__main__":
    from sympy import pretty_print
    from schubmult.abc import x,y
    from schubmult.rings.polynomial_algebra import ForestPolyBasis, PolynomialAlgebra
    ForestPoly = PolynomialAlgebra(ForestPolyBasis(x))
    ForestPolyY = PolynomialAlgebra(ForestPolyBasis(y))
    poly = double_forest_polynomial([2, 0, 1, 1], lambda i: x[i], lambda j: y[j])
    first_forest = ForestPoly.from_expr(poly,length=4)

    result = 0
    pretty_print(poly)
    for comp, coeff in first_forest.items():
        result += ForestPolyY.from_expr(coeff,length=4) @ ForestPolyY(*comp)
    pretty_print(result)
