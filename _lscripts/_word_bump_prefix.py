"""Word-level tests of the Little bump chain.
(A) Is the chain of bumped positions monotone in position?
(B) For w = u v reduced and a bump at any position q, is the restriction of the chain to u empty or equal to the
    Little bump of u alone at the entry position (largest chain position inside u)?"""
import random
from itertools import permutations
from schubmult import Permutation
from schubmult.utils.perm_utils import find_reduced_fail, is_reduced
from collections import Counter

def ltb_chain(word, index):
    word = [*word]
    chain = []
    while True:
        if word[index] == 1:
            return None, None
        word[index] -= 1
        chain.append(index)
        if is_reduced(word):
            return tuple(word), chain
        index = find_reduced_fail(word, index)
        if index is None:
            return None, None

def reduced_words(w):
    # all reduced words of w via recursion on descents
    if w.inv == 0:
        yield ()
        return
    for d in w.descents():
        # right descent d (0-indexed?) -> letter d+1
        sw = w * Permutation.ref_product(d + 1) if hasattr(Permutation, "ref_product") else None
        for rw in reduced_words(sw):
            yield (*rw, d + 1)

stats = Counter()
bad = []
random.seed(1)
for N in range(3, 7):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv == 0 or w.inv > 7:
            continue
        words = list(reduced_words(w))
        if len(words) > 40:
            words = random.sample(words, 40)
        for word in words:
            L = len(word)
            for q in range(L):
                new, chain = ltb_chain(word, q)
                if new is None:
                    continue
                # (A) monotone?
                mono_dec = all(chain[k] > chain[k + 1] for k in range(len(chain) - 1))
                mono_inc = all(chain[k] < chain[k + 1] for k in range(len(chain) - 1))
                stats[("chain", "dec" if mono_dec else "inc" if mono_inc else "nonmono")] += 1
                # (B) prefix restriction for every cut k
                for k in range(1, L):
                    u = word[:k]
                    u_new = new[:k]
                    inside = [c for c in chain if c < k]
                    if not inside:
                        stats[("B", "unchanged", u_new == u)] += 1
                        continue
                    entry = max(inside)
                    bu, _ = ltb_chain(u, entry)
                    ok = bu == u_new
                    stats[("B", "entry-bump", ok)] += 1
                    if not ok and len(bad) < 5:
                        bad.append((word, q, chain, k, u_new, bu))
print(stats)
for b in bad:
    print(b)
