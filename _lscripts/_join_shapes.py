"""Local confluence shapes for pairs (beta, L) where L = canonical bump (last descent (h,h+1) or Little (r,s)).
Classify the minimal join: L(beta T) == L(T); L(beta T) == beta'(L T); L^2(beta T) == beta' L T; L(beta T)== beta'beta'' L T; etc."""
import random, sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import rroot, ltb_down, reduced_words, maxd, valid_rc
from collections import Counter

def perm_arr(word, n=12):
    w = Permutation.ref_product(*word) if word else Permutation([])
    arr = list(w); arr += list(range(len(arr) + 1, n + 1))
    return arr[:n]

def high_bumps(word, p):
    out = {}
    for q in range(len(word)):
        a, b = rroot(word, q)
        if a > p:
            nw = ltb_down(word, q)
            if nw is not None:
                out[(a, b)] = nw
    return out

def canon(word, p, mode):
    """Canonical bump: mode 'desc' -> root (h,h+1); mode 'little' -> lex largest inversion (r,s)."""
    arr = perm_arr(word)
    h = maxd(Permutation.ref_product(*word))
    if h <= p:
        return None
    if mode == "desc":
        root = (h, h + 1)
    else:
        s = max(j for j in range(h + 1, len(arr) + 1) if arr[j - 1] < arr[h - 1])
        root = (h, s)
    for q in range(len(word)):
        if rroot(word, q) == root:
            return ltb_down(word, q)
    return None

def canon_iter(word, p, mode, k):
    for _ in range(k):
        if word is None:
            return None
        word = canon(word, p, mode)
    return word

random.seed(0)
for mode in ("desc", "little"):
    stats = Counter()
    for N in range(3, 8):
        for arr in permutations(range(1, N + 1)):
            w = Permutation(list(arr))
            if w.inv == 0 or w.inv > 8:
                continue
            words = list(reduced_words(w))
            if len(words) > 12:
                words = random.sample(words, 12)
            for word in words:
                for p in range(1, maxd(w)):
                    L = canon(word, p, mode)
                    if L is None:
                        continue
                    for key, bT in high_bumps(word, p).items():
                        if bT == L:
                            continue
                        LbT = canon(bT, p, mode)
                        if LbT is None:
                            stats[(mode, "beta T already normal")] += 1
                            continue
                        if LbT == L:
                            stats[(mode, "L(bT)=L(T)")] += 1
                            continue
                        hbL = set(high_bumps(L, p).values())
                        if LbT in hbL:
                            stats[(mode, "L(bT)=b'(LT)")] += 1
                            continue
                        # deeper
                        L2bT = canon(LbT, p, mode)
                        if L2bT is not None and (L2bT == L or L2bT in hbL):
                            stats[(mode, "L^2(bT) in {LT, b'LT}")] += 1
                            continue
                        hb2L = set()
                        for x in hbL:
                            hb2L |= set(high_bumps(x, p).values())
                        if LbT in hb2L:
                            stats[(mode, "L(bT)=b'b''(LT)")] += 1
                            continue
                        if L2bT is not None and L2bT in hb2L:
                            stats[(mode, "L^2(bT)=b'b''(LT)")] += 1
                            continue
                        stats[(mode, "OTHER")] += 1
    print(stats)
