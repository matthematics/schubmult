"""Unsigned marginal counts for ordinary Schubert times Schur (not a fast multiplication kernel).

coefficient(u, shape, k, w) computes one stable coefficient. product(u, shape, k)
computes the stable product, not its truncation to the input symmetric group.
The construction uses Pawlowski's positroid diagram theorem (arXiv:1612.00097,
Theorem 1.6 and Proposition 7.10). See the companion note schubert-schur-margins.tex.
"""

from __future__ import annotations

from collections import defaultdict
from itertools import pairwise
from math import comb, factorial
from typing import TYPE_CHECKING

from schubmult.combinatorics.permutation import Permutation

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence

Cell = tuple[int, int]
Word = tuple[int, ...]
Bin = frozenset[int]


def _shape(shape: Sequence[int], k: int) -> Word:
    if k < 1 or any(a < 0 for a in shape) or any(a < b for a, b in pairwise(shape)):
        raise ValueError("Expected k >= 1 and a weakly decreasing nonnegative partition")
    positive = tuple(a for a in shape if a)
    if len(positive) > k:
        raise ValueError("The partition has more than k parts")
    return positive + (0,) * (k - len(positive))


def _classes(keys: Iterable[int]) -> list[list[int]]:
    groups: dict[int, list[int]] = defaultdict(list)
    for i, key in enumerate(keys):
        groups[key].append(i)
    return list(groups.values())


def _bins(cells: Sequence[Cell], mu: Word) -> tuple[int, list[Bin], list[Bin]]:
    """Repeated letters in a column have zero coordinate, so omit them."""
    rows = _classes(r for r, _ in cells)
    cols = _classes(c for _, c in cells)
    words: list[Word] = []
    counts = list(mu)
    used: dict[int, set[int]] = defaultdict(set)
    prefix: list[int] = []

    def visit(i: int) -> None:
        if i == len(cells):
            words.append(tuple(prefix))
            return
        col = cells[i][1]
        for letter in range(len(mu)):
            if counts[letter] and letter not in used[col]:
                prefix.append(letter)
                counts[letter] -= 1
                used[col].add(letter)
                visit(i + 1)
                used[col].remove(letter)
                counts[letter] += 1
                prefix.pop()

    visit(0)
    index = {w: i for i, w in enumerate(words)}
    fixed: set[Bin] = set()
    raising: dict[tuple[int, Word], set[int]] = defaultdict(set)
    by_row: dict[tuple[Word, ...], set[int]] = defaultdict(set)
    for i, w in enumerate(words):
        for col in cols:
            for a, b in pairwise(col):
                other = list(w)
                other[a], other[b] = other[b], other[a]
                fixed.add(frozenset((i, index[tuple(other)])))
        for pos, letter in enumerate(w):
            if letter:
                raised = (*w[:pos], letter - 1, *w[pos + 1 :])
                raising[letter, raised].add(i)
        row_type = tuple(tuple(sorted(w[j] for j in row)) for row in rows)
        by_row[row_type].add(i)
    fixed.update(frozenset(b) for b in raising.values())
    return len(words), sorted(fixed, key=lambda b: (len(b), sorted(b))), [frozenset(by_row[key]) for key in sorted(by_row)]


def _marked(n: int, fixed: list[Bin], readouts: list[Bin]) -> list[int]:
    """Balance decks one bin at a time, using only unions and copies (in compressed form)."""
    bins = fixed + readouts
    empty = (0,) * len(bins)
    decks = [(tuple(int(i in h) for h in bins), empty) for i in range(n)]
    # Optional exact compression: all relevant 0/1 minors have absolute value <= n!.
    # Use natural integers if none of these known Mersenne primes certifies the computation.
    bound = factorial(n)
    primes = ((1 << e) - 1 for e in (31, 61, 127, 521, 2203, 19937))
    prime = next((m for m in primes if m > bound), None)
    marked = []
    for h in range(len(bins)):
        pos = next((i for i, (left, right) in enumerate(decks) if left[0] != right[0]), None)
        if pos is None:
            decks = [(p[1:], q[1:]) for p, q in decks]
            continue
        if h >= len(fixed):
            marked.append(h - len(fixed))
        p, q = decks.pop(pos)
        a, b = p[0], q[0]
        next_decks = []
        for r, s in decks:
            c, d = r[0], s[0]
            if c == d:
                next_decks.append((r[1:], s[1:]))
                continue
            left = tuple(c * x + d * y + b * z + a * t for x, y, z, t in zip(p[1:], q[1:], r[1:], s[1:]))
            right = tuple(c * y + d * x + b * t + a * z for x, y, z, t in zip(p[1:], q[1:], r[1:], s[1:]))
            if prime is not None:
                left = tuple(x % prime for x in left)
                right = tuple(x % prime for x in right)
            next_decks.append((left, right))
        decks = next_decks
        if not decks:
            break
    return marked


def hook_count(cells: Sequence[Cell], leg: int) -> int:
    """Count subsets of non-root components of a forest for a complementary hook."""
    vertices: dict[tuple[str, int], tuple[str, int]] = {}

    def root(v: tuple[str, int]) -> tuple[str, int]:
        vertices.setdefault(v, v)
        while vertices[v] != v:
            vertices[v] = vertices[vertices[v]]
            v = vertices[v]
        return v

    for r, c in cells:
        a, b = root(("r", r)), root(("c", c))
        if a == b:
            return 0
        vertices[a] = b
    extra = len(cells) - len({c for _, c in cells})
    components = len({root(v) for v in vertices})
    return comb(components - 1, leg - extra) if extra <= leg < extra + components else 0


def diagram_count(cells: Sequence[Cell], mu: Sequence[int]) -> int:
    """Multiplicity for a diagram and a partition of its size."""
    if len(set(cells)) != len(cells):
        raise ValueError("Repeated diagram cell")
    mu = tuple(a for a in _shape(mu, max(1, len(mu))) if a)
    if sum(mu) != len(cells):
        raise ValueError("Diagram and partition sizes differ")
    cells = tuple(sorted(cells))
    while mu:
        if len(mu) == 1 or mu[1] == 1:
            return hook_count(cells, len(mu) - 1)
        columns = _classes(c for _, c in cells)
        if any(len(col) > len(mu) for col in columns):
            return 0
        full = {cells[col[0]][1] for col in columns if len(col) == len(mu)}
        if not full:
            return len(_marked(*_bins(cells, mu)))
        if len(full) > mu[-1]:
            return 0
        cells = tuple((r, c) for r, c in cells if c not in full)
        mu = tuple(a - len(full) for a in mu if a > len(full))
    return 1


def coefficient(u: Sequence[int], shape: Sequence[int], k: int, w: Sequence[int]) -> int:
    """Stable ordinary coefficient of S_w in S_u s_shape(x_1,...,x_k)."""
    lam = _shape(shape, k)
    pu, pw = Permutation(u), Permutation(w)
    if pw.inv != pu.inv + sum(lam):
        return 0
    if not lam[0]:
        return int(pu == pw)
    n = max(len(pu), len(pw), k + lam[0])
    x, y = tuple(pu[i] for i in range(n)), tuple(pw[i] for i in range(n))
    if any(x[i] > y[i] for i in range(k)) or any(x[i] < y[i] for i in range(k, n)):
        return 0
    if any(x[i] < x[j] and y[i] > y[j] for block in (range(k), range(k, n)) for i in block for j in block if i < j):
        return 0
    f = [0] * n
    for i in range(n):
        f[y[i] - 1] = x[i] + n * (i < k)
    cells = [(i + 1, 1 + (f[j % n] - 1) % n) for i in range(n) for j in range(i + 1, i + n) if f[i] > f[j % n] + n * (j // n)]
    return diagram_count(cells, tuple(n - k - a for a in reversed(lam)))


def product(u: Sequence[int], shape: Sequence[int], k: int) -> dict[Permutation, int]:
    """Use the support of iterated column Pieri as a finite list of candidate endpoints."""
    from schubmult.utils.schub_lib import elem_sym_perms

    lam = _shape(shape, k)
    start = Permutation(u)
    endpoints = {start}
    for j in range(lam[0]):
        height = sum(a > j for a in lam)
        endpoints = {w for p in endpoints for w, steps in elem_sym_perms(p, height, k) if steps == height}
    return {w: c for w in endpoints if (c := coefficient(tuple(start), lam, k, tuple(w)))}
