"""Which invariants do p-high bumps preserve? Test w|_{[p]} (values at positions <= p) and the set of those values.
Also whether NF_p(word) is determined by (Q, w|_{[p]}) or by (Q, w restricted to positions <= p with flattening)."""
import random
from itertools import permutations
from schubmult import Permutation
from schubmult.utils.perm_utils import find_reduced_fail, is_reduced
from collections import Counter, defaultdict

def rroot(word, q):
    a, b = word[q], word[q] + 1
    for x in word[q + 1:]:
        a = a + 1 if a == x else a - 1 if a == x + 1 else a
        b = b + 1 if b == x else b - 1 if b == x + 1 else b
    return (min(a, b), max(a, b))

def ltb_down(word, index):
    word = [*word]
    while True:
        if word[index] == 1:
            return None
        word[index] -= 1
        if is_reduced(word):
            return tuple(word)
        index = find_reduced_fail(word, index)
        if index is None:
            return None

def high_bumps(word, p):
    out = set()
    for q in range(len(word)):
        a, b = rroot(word, q)
        if a > p:
            nw = ltb_down(word, q)
            if nw is not None:
                out.add(nw)
    return out

memo = {}
def NF(word, p):
    key = (word, p)
    if key in memo:
        return memo[key]
    nb = high_bumps(word, p)
    res = word if not nb else NF(next(iter(nb)), p)
    memo[key] = res
    return res

def perm_of(word, n):
    w = Permutation.ref_product(*word) if word else Permutation([])
    arr = list(w)
    arr += list(range(len(arr) + 1, n + 1))
    return tuple(arr[:n])

def eg_Q(word):
    P, Q = [], []
    for t, x in enumerate(word, 1):
        r = 0
        while True:
            if r == len(P):
                P.append([x]); Q.append([t]); break
            row = P[r]
            bigger = [k for k, y in enumerate(row) if y > x]
            if not bigger:
                row.append(x); Q[r].append(t); break
            k = bigger[0]; y = row[k]
            if y == x + 1 and k > 0 and row[k - 1] == x:
                x = y
            else:
                row[k] = x; x = y
            r += 1
    return tuple(tuple(r) for r in Q)

def reduced_words(w):
    if w.inv == 0:
        yield ()
        return
    for d in w.descents():
        sw = w * Permutation.ref_product(d + 1)
        for rw in reduced_words(sw):
            yield (*rw, d + 1)

random.seed(0)
pres = Counter()
tbl = defaultdict(set)
n = 9
for N in range(3, 7):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv == 0 or w.inv > 7:
            continue
        words = list(reduced_words(w))
        if len(words) > 30:
            words = random.sample(words, 30)
        for word in words:
            for p in range(1, len(w.trimcode)):
                wp = perm_of(word, n)
                for nw in high_bumps(word, p):
                    wq = perm_of(nw, n)
                    pres[("values<=p same", wp[:p] == wq[:p])] += 1
                    pres[("set same", set(wp[:p]) == set(wq[:p]))] += 1
                nf = NF(word, p)
                tbl[(p, eg_Q(word), wp[:p])].add(nf)
print(pres)
amb = {k: v for k, v in tbl.items() if len(v) > 1}
print("keys", len(tbl), "ambiguous (Q, values at positions<=p):", len(amb))
for k, v in list(amb.items())[:4]:
    print(k, sorted(v))
