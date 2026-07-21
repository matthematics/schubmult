"""Independent test of the LITERAL paper statement (Theorem 1.2, quantum Pieri) of
"A Molev-Sagan type formula for quantum double Schubert polynomials".

Unlike verify_quantum_pieri.py (which only checks the ring code is self-consistent),
this implements Definition 1.1/1.2 combinatorics DIRECTLY from the paper text:

    S^q_u(x;y) E^q_{p,k}(x;z)
        = sum_{u ->_k w} q_k(u,w) * E_{p-n_k(u,w), k-n_k(u,w)}(y_{P_k(u,w)}; z) * S^q_w(x;y)

Definitions implemented from scratch:
  u ->_k w : exists reflections t_{a_i b_i} with a_i<=k<b_i, a_i distinct,
             b_i nonincreasing, each step d ell in {+1, -2(b_i-a_i)+1}, w = u * prod.
  q_k(u,w): along the sequence, multiply by q_{a b} = q_a q_{a+1} ... q_{b-1}
             exactly on the length-DROP steps.
  P_k(u,w) = { u(i) : i<=k and u(i)=w(i) }        (fixed values in first k positions)
  n_k(u,w) = #{ i<=k : u(i) != w(i) }              ( = k - |P_k| )

RHS uses the NON-quantum factorial elementary E = FactorialElemSym in vars y_{P_k}, z.

We also verify q_k(u,w) is well-defined (all reflection sequences u->_k w agree).

Run: conda activate schubmult_312 && python _lscripts/verify_paper_pieri.py 4
"""

import sys

from schubmult import *  # noqa: F401,F403  (init package to avoid circular import)
from schubmult.abc import x, z
from schubmult.combinatorics.permutation import Permutation
from schubmult.rings.schubert import QDSx
from schubmult.symbolic import S, expand, expand_func, prod
from schubmult.symbolic.poly.variables import GeneratingSet
from schubmult.symbolic.symmetric_polynomials import FactorialElemSym

y = GeneratingSet("y")
z_gs = GeneratingSet("z")
q_gs = GeneratingSet("q")


def length(perm):
    """Number of inversions of a Permutation (one-line, 1-based values)."""
    w = list(perm)
    n = len(w)
    inv = 0
    for i in range(n):
        for j in range(i + 1, n):
            if w[i] > w[j]:
                inv += 1
    return inv


def apply_t(perm, a, b):
    """Right-multiply by transposition t_{a,b}: swap the entries in positions a,b (1-based)."""
    return perm.swap(a - 1, b - 1)


def q_ab(a, b):
    """q_{a,b} = q_a q_{a+1} ... q_{b-1}."""
    return prod([q_gs[s] for s in range(a, b)])


def enumerate_pieri(u, k, N):
    """Return dict {w: qweight} for all w with u ->_k w, plus a consistency flag.

    BFS over reflection sequences per Definition 1.2. State carries used a-set,
    the max allowed b (nonincreasing constraint), and the accumulated q-weight.
    """
    results = {}  # w (tuple) -> set of qweights seen
    u = Permutation(u)
    lu = length(u)
    # state: (perm, used_a frozenset, last_b, qweight, curr_len)
    start = (u, frozenset(), N + 1, S.One, lu)
    stack = [start]
    # record start (p=0, w=u)
    results.setdefault(tuple(u), set()).add(S.One)
    seen_states = set()
    while stack:
        perm, used_a, last_b, qw, clen = stack.pop()
        state_key = (tuple(perm), used_a, last_b)
        if state_key in seen_states:
            continue
        seen_states.add(state_key)
        for a in range(1, k + 1):
            if a in used_a:
                continue
            for b in range(k + 1, N + 1):
                if b > last_b:
                    continue
                nperm = apply_t(perm, a, b)
                nlen = length(nperm)
                d = nlen - clen
                if d == 1:
                    nqw = qw
                elif d == -2 * (b - a) + 1:
                    nqw = qw * q_ab(a, b)
                else:
                    continue
                results.setdefault(tuple(nperm), set()).add(nqw)
                stack.append((nperm, used_a | {a}, b, nqw, nlen))
    consistent = all(len(v) == 1 for v in results.values())
    weights = {w: next(iter(v)) for w, v in results.items()}
    return weights, consistent


def paper_rhs(u, p, k, N, ring):
    u = Permutation(u)
    weights, consistent = enumerate_pieri(u, k, N)
    rhs = ring.zero
    for w_tuple, qw in weights.items():
        w = Permutation(w_tuple)
        # P_k and n_k
        Pk = [u[i] for i in range(k) if i < len(u) and u[i] == (w[i] if i < len(w) else i + 1)]
        # handle positions beyond len via identity value i+1
        Pk = []
        for i in range(k):
            uv = u[i] if i < len(u) else i + 1
            wv = w[i] if i < len(w) else i + 1
            if uv == wv:
                Pk.append(uv)
        nk = k - len(Pk)
        deg = p - nk
        nv = k - nk
        if deg < 0 or deg > nv:
            continue  # elementary is zero
        yvars = [y[v] for v in sorted(Pk)]
        ncoeff = nv + 1 - deg  # number of z coeff vars needed
        zvars = [z[i] for i in range(1, ncoeff + 1)]
        Eterm = FactorialElemSym(deg, nv, yvars, zvars) if nv > 0 else (S.One if deg == 0 else S.Zero)
        rhs += (qw * expand_func(Eterm)) * ring(w)
    return rhs, consistent


def cycle_cpk(p, k):
    return Permutation.ref_product(*range(k - p + 1, k + 1))


def elems_equal(a, b):
    diff = a - b
    for v in diff.values():
        if expand(v) != 0:
            return False
    return True


def run(n):
    N = n + n  # generous room for support beyond n
    perms = list(Permutation.all_permutations(n))
    total = 0
    mismatches = []
    inconsistent = []
    for u in perms:
        Su = QDSx(u)
        ring = Su.ring
        for k in range(1, n):
            for p in range(1, k + 1):
                c = cycle_cpk(p, k)
                lhs = Su * QDSx(c, z_gs)
                rhs, consistent = paper_rhs(u, p, k, N, ring)
                total += 1
                if not consistent:
                    inconsistent.append((tuple(u), p, k))
                if not elems_equal(lhs, rhs):
                    mismatches.append((tuple(u), p, k))
    print(f"n={n}: tested {total} (u,p,k); mismatches={len(mismatches)}; q_k-inconsistent={len(inconsistent)}")
    for m in mismatches[:15]:
        print("  MISMATCH:", m)
    for m in inconsistent[:15]:
        print("  Q_K INCONSISTENT:", m)
    return len(mismatches) == 0 and len(inconsistent) == 0


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    ok = run(n)
    sys.exit(0 if ok else 1)
