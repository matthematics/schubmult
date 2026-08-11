"""Independent check of the counterexample to  N(u, mu_v) <= n^k N(u,v).

Compares the union-of-supports count with a genuine two-alphabet double product
S_u(x;y) * S_v(x;z), counting Schubert-basis terms with nonzero coefficient.
"""
from schubmult import Permutation
from schubmult.rings.schubert.schubert_ring import Sx
from schubmult.rings.schubert.double_schubert_ring import DSx
from schubmult.symbolic import sympify


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
    return [Permutation(list(t)) for t in seen]


def N_union(u, w, n):
    acc = set()
    for wp in left_weak_downset(w, n):
        acc |= set((Sx(u) * Sx(wp)).keys())
    return len(acc)


def N_double(u, w):
    prod = DSx(u, "y") * DSx(w, "z")
    cnt = 0
    for _, coeff in prod.items():
        if sympify(coeff) != 0:
            cnt += 1
    return cnt


n = 5
u = Permutation([2, 1, 3, 5, 4])
v = Permutation([1, 2, 5, 3, 4])
mu = ~((~v).minimal_dominant_above())
k = len(v.descents())
print("u=", list(u), "v=", list(v), "mu=", list(mu), "k=", k, "n^k=", n**k)
print("N_union(u,v)   =", N_union(u, v, n), "  N_double(u,v)   =", N_double(u, v))
print("N_union(u,mu)  =", N_union(u, mu, n), "  N_double(u,mu)  =", N_double(u, mu))
print("interval |[e,v]_L| =", len(left_weak_downset(v, n)),
      " |[e,mu]_L| =", len(left_weak_downset(mu, n)))
