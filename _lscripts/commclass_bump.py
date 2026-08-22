"""Test the key lemma: bumping is a well-defined operation on commutation classes.

Claim: if w1, w2 are reduced words in the same commutation class, and
t_1..t_a, t'_1..t'_b are the bumping sequences produced by iterating the Little
bump at the last descent, then a = b and for every k the words
w1 up_{t_1..t_k} and w2 up_{t'_1..t'_k} are reduced words for the same permutation
(indeed, empirically, lie in the same commutation class).
"""

import itertools
import sys
from collections import Counter, deque

from schubmult import Permutation


def perm_of(word):
    p = Permutation([])
    for a in word:
        p = p.swap(a - 1, a)
    return p


def is_reduced(word):
    return perm_of(word).inv == len(word)


def commutation_class(word):
    """All reduced words obtained from `word` by adjacent commutations |a-b|>=2."""
    seen = {tuple(word)}
    queue = deque([tuple(word)])
    while queue:
        w = queue.popleft()
        for i in range(len(w) - 1):
            if abs(w[i] - w[i + 1]) >= 2:
                v = list(w)
                v[i], v[i + 1] = v[i + 1], v[i]
                v = tuple(v)
                if v not in seen:
                    seen.add(v)
                    queue.append(v)
    return seen


def little_bump(word, pos):
    """Little bump of `word` beginning at position `pos` (0-indexed)."""
    w = list(word)
    i = pos
    while True:
        if w[i] > 1:
            w[i] -= 1
        else:
            w = [x + 1 if j != i else x for j, x in enumerate(w)]
        if is_reduced(w):
            return tuple(w)
        # unique j != i whose deletion restores reducedness
        cand = [j for j in range(len(w)) if j != i and is_reduced(w[:j] + w[j + 1:])]
        if len(cand) != 1:
            raise RuntimeError(f"non-unique defect in {w} at {i}: {cand}")
        i = cand[0]


def last_descent_pos(word):
    """Position of the crossing at the last descent of the permutation."""
    p = perm_of(word)
    if p.inv == 0:
        return None
    d = max(p.descents())
    # the last letter of the word equal to the last descent index is bumped
    # (Little's algorithm: swap location of the last inversion)
    m = max(word)
    idx = [i for i, a in enumerate(word) if a == m]
    return idx[-1] if idx else None


def bump_chain(word):
    """Iterate the bump at the last descent until max letter drops; record words."""
    chain = [tuple(word)]
    w = tuple(word)
    guard = 0
    while True:
        guard += 1
        if guard > 50:
            raise RuntimeError("no termination")
        p = perm_of(w)
        if p.inv == 0:
            return chain
        pos = last_descent_pos(w)
        if pos is None:
            return chain
        target = max(w)
        w = little_bump(w, pos)
        chain.append(w)
        if max(w) < target:
            return chain


def all_reduced_words(perm):
    """All reduced words of `perm`, by peeling descents."""
    if perm.inv == 0:
        return [()]
    out = []
    for d in perm.descents():
        sub = perm.swap(d, d + 1)
        for w in all_reduced_words(sub):
            out.append((*w, d + 1))
    return list(dict.fromkeys(out))


def main(n):
    stats = Counter()
    failures = []
    seen_classes = set()

    for perm in Permutation.all_permutations(n):
        if perm.inv == 0:
            continue
        for word in perm.get_reduced_words():
            cls = frozenset(commutation_class(word))
            if cls in seen_classes:
                continue
            seen_classes.add(cls)
            reps = sorted(cls)
            if len(reps) < 2:
                continue
            chains = [bump_chain(w) for w in reps]
            lengths = {len(c) for c in chains}
            stats["classes"] += 1
            if len(lengths) != 1:
                stats["LENGTH_MISMATCH"] += 1
                if len(failures) < 3:
                    failures.append(("length", reps, [len(c) for c in chains]))
                continue
            stats["same_length"] += 1
            ok_perm = True
            ok_class = True
            for k in range(len(chains[0])):
                perms = {perm_of(c[k]) for c in chains}
                if len(perms) != 1:
                    ok_perm = False
                classes = {frozenset(commutation_class(c[k])) for c in chains}
                if len(classes) != 1:
                    ok_class = False
            stats["same_perm_each_step" if ok_perm else "PERM_MISMATCH"] += 1
            stats["same_class_each_step" if ok_class else "class_differs"] += 1
            if not ok_perm and len(failures) < 3:
                failures.append(("perm", reps, None))

    print(f"n = {n}")
    for k in sorted(stats):
        print(f"  {k}: {stats[k]}")
    for f in failures:
        print("  failure:", f)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
