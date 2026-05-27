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
    """Weight map in the vine model."""
    i = letter["syllable"]
    j = letter["value"]
    if letter["barred"]:
        return t_gen(j) - t_gen(i)
    return x_gen(i) - t_gen(j)


def double_forest_polynomial(code, x_gen, t_gen, n=None):
    """Compute P_F(x;t) for indexed forest F with code c(F)=code."""
    code = list(code)
    while code and code[-1] == 0:
        code.pop()
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
