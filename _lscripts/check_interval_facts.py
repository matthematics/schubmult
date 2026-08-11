"""Check the interval-level facts underpinning the salvageable results:
  (A)  Lambda(mu_v) = |[e,mu_v]_L| <= n^k |[e,v]_L|        (k=#descents v)
  (B)  |[e,v]_L| <= N(u,v)  for all u                      (so (A) => Lambda<=n^k N)
"""
import sys
from schubmult import Permutation
from schubmult.rings.schubert.schubert_ring import Sx


def left_weak_downset(perm, n):
    start = tuple(perm[i] for i in range(n))
    seen = {start}
    stack = [start]
    while stack:
        a = stack.pop()
        pos = [0] * (n + 2)
        for idx, val in enumerate(a):
            pos[val] = idx
        for i in range(1, n):
            if pos[i] > pos[i + 1]:
                b = list(a)
                b[pos[i]] = i + 1
                b[pos[i + 1]] = i
                b = tuple(b)
                if b not in seen:
                    seen.add(b)
                    stack.append(b)
    return seen


def check_A(n):
    worst = 0.0
    info = None
    for v in Permutation.all_permutations(n):
        k = len(v.descents())
        mu = ~((~v).minimal_dominant_above())
        lv = len(left_weak_downset(v, n))
        lmu = len(left_weak_downset(mu, n))
        r = lmu / (n ** k * lv)
        if r > worst:
            worst = r
            info = (tuple(v), tuple(mu), k, lv, lmu, n ** k)
    print(f"(A) n={n}: max |[e,mu]_L|/(n^k |[e,v]_L|) = {worst:.4f}  {info}")


def check_B(n):
    perms = Permutation.all_permutations(n)
    worst = 0.0
    info = None
    for v in perms:
        v_ds = [Permutation(list(t)) for t in left_weak_downset(v, n)]
        lv = len(v_ds)
        for u in perms:
            acc = set()
            for wp in v_ds:
                acc |= set((Sx(u) * Sx(wp)).keys())
            Nv = len(acc)
            r = lv / Nv
            if r > worst:
                worst = r
                info = (tuple(u), tuple(v), lv, Nv)
    print(f"(B) n={n}: max |[e,v]_L|/N(u,v) = {worst:.4f}  {info}  (want <=1)")


if __name__ == "__main__":
    for n in range(2, 8):
        check_A(n)
    for n in (4, 5):
        check_B(n)
