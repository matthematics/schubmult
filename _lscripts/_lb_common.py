"""Shared helpers for Little bump experiments (non-degenerate chains only)."""
from schubmult import Permutation
from schubmult.utils.perm_utils import is_reduced

def rroots(word):
    """Right-transported roots (positions convention) of every letter, as ordered pairs (may be negative)."""
    n = len(word)
    out = [None] * n
    for q in range(n):
        a, b = word[q], word[q] + 1
        for x in word[q + 1:]:
            a = a + 1 if a == x else a - 1 if a == x + 1 else a
            b = b + 1 if b == x else b - 1 if b == x + 1 else b
        out[q] = (a, b)
    return out

def rroot(word, q):
    a, b = rroots(word)[q]
    return (min(a, b), max(a, b))

def ltb_down(word, index, return_chain=False):
    """Downward Little bump at ``index``; None if degenerate (letter 1, or non-unique partner)."""
    word = [*word]
    chain = []
    while True:
        if word[index] == 1:
            return (None, None) if return_chain else None
        word[index] -= 1
        chain.append(index)
        if is_reduced(word):
            return (tuple(word), chain) if return_chain else tuple(word)
        rts = rroots(word)
        target = set(rts[index])
        partners = [i for i in range(len(word)) if i != index and set(rts[i]) == target]
        if len(partners) != 1:
            return (None, None) if return_chain else None
        index = partners[0]

def valid_rc(word, seq):
    for k in range(len(word)):
        if word[k] < seq[k]:
            return False
        if k + 1 < len(word) and seq[k] == seq[k + 1] and not word[k] > word[k + 1]:
            return False
    return True

def reduced_words(w):
    if w.inv == 0:
        yield ()
        return
    for d in w.descents():
        sw = w * Permutation.ref_product(d + 1)
        for rw in reduced_words(sw):
            yield (*rw, d + 1)

def maxd(w):
    return len(w.trimcode)
