"""Stress-test the per-case structural claims used in the Proof of Theorem 1.2.

The inductive proof extracts the coefficient of S^q_w from the recursion (Lemma 3.12)
        S_u E^q_{p,k} = (x_k - z_{k-p+1}) S_u E^q_{p-1,k-1} + S_u E^q_{p,k-1} + q_{k-1} S_u E^q_{p-2,k-2}
by applying the single-variable quantum Monk (Lemma 3.2) at position i=k to the first term.
For a fixed target w, the Monk step reaches w from a source v=w t_{k,pp} (pp != k), where the
transition length type is RAISE (dl=+1) or DROP (dl=1-2|pp-k|); the source contributes iff
u ->_{k-1} v.  We group these by pp>k (UP, sign +) and pp<k (DOWN, sign -):

    UP(u,w,k)   = { b>k : v=w t_{kb}, dl(w<-v) in {+1, 1-2(b-k)}, u ->_{k-1} v }
    DOWN(u,w,k) = { a<k : v=w t_{ak}, dl(w<-v) in {+1, 1-2(k-a)}, u ->_{k-1} v }

This script independently implements ->_m, q_m, P_m, n_m (Definition 1.1) via BFS and checks
the exact claims the proof requires of Lemmas 3.8-3.11, split by the four cases

    (1) u ->_k w  and  u ->_{k-1} w         (Lemma 3.9,  "pieritogether")
    (2) u ->_k w  and  u -/->_{k-1} w        (Lemma 3.8,  "pieriknotkn1")
    (3) u -/-> _k w and  u ->_{k-1} w        (Lemma 3.10, "pierikn1notk")
    (4) u -/-> _k w and  u -/->_{k-1} w      (Lemma 3.11, "pieribothnot")

For each case it reports violations of the precise structural claim, so we learn the exact
correct statements to prove.

Run: conda activate schubmult_312 && python _lscripts/stress_test_pieri_cases.py 5
"""

import sys
from collections import defaultdict

from schubmult import *  # noqa: F401,F403
from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import S, prod
from schubmult.symbolic.poly.variables import GeneratingSet

q_gs = GeneratingSet("q")


def q_ab(a, b):
    return prod([q_gs[s] for s in range(a, b)])


def val(perm, pos):
    """1-indexed value perm(pos), using identity beyond stored support."""
    return perm[pos - 1]


def P_set(u, w, m):
    return frozenset(val(u, i) for i in range(1, m + 1) if val(u, i) == val(w, i))


def n_val(u, w, m):
    return sum(1 for i in range(1, m + 1) if val(u, i) != val(w, i))


def enumerate_pieri(u, k, N):
    """Return {w(Permutation): q_k(u,w)} for all w with u ->_k w (Definition 1.1)."""
    u = Permutation(u)
    results = {u: S.One}
    start = (u, frozenset(), N + 1, S.One, u.inv)
    stack = [start]
    seen = set()
    while stack:
        perm, used_a, last_b, qw, clen = stack.pop()
        key = (perm, used_a, last_b)
        if key in seen:
            continue
        seen.add(key)
        for a in range(1, k + 1):
            if a in used_a:
                continue
            for b in range(k + 1, N + 1):
                if b > last_b:
                    continue
                nperm = perm.swap(a - 1, b - 1)
                d = nperm.inv - clen
                if d == 1:
                    nqw = qw
                elif d == -2 * (b - a) + 1:
                    nqw = qw * q_ab(a, b)
                else:
                    continue
                # record smallest (consistent) weight; consistency verified elsewhere
                results.setdefault(nperm, nqw)
                stack.append((nperm, used_a | {a}, b, nqw, nperm.inv))
    return results


def transition_type(w, v, lo, hi):
    """dl = inv(w) - inv(v); return 'raise' if +1, 'drop' if 1-2(hi-lo), else None."""
    dl = w.inv - v.inv
    if dl == 1:
        return "raise"
    if dl == 1 - 2 * (hi - lo):
        return "drop"
    return None


def up_down(u, w, k, R_km1, N):
    """Return (UP, DOWN) lists of (partner, type, v) with u ->_{k-1} v."""
    up, down = [], []
    for b in range(k + 1, N + 1):
        v = w.swap(k - 1, b - 1)
        if v not in R_km1:
            continue
        t = transition_type(w, v, k, b)
        if t is not None:
            up.append((b, t, v))
    for a in range(1, k):
        v = w.swap(a - 1, k - 1)
        if v not in R_km1:
            continue
        t = transition_type(w, v, a, k)
        if t is not None:
            down.append((a, t, v))
    return up, down


def run(n):
    N = 2 * n
    perms = list(Permutation.all_permutations(n))
    # counters
    counts = defaultdict(int)
    viol = defaultdict(list)

    for u in perms:
        for k in range(2, n):  # k>1; k=1 is the base case
            R_k = enumerate_pieri(u, k, N)
            R_km1 = enumerate_pieri(u, k - 1, N)
            R_km2 = enumerate_pieri(u, k - 2, N) if k - 2 >= 0 else {}

            # candidate targets: anything with a nonzero contribution on either side
            cands = set(R_k) | set(R_km1) | set(R_km2)
            for v in list(R_km1):
                for pp in range(1, N + 1):
                    if pp == k:
                        continue
                    cands.add(v.swap(k - 1, pp - 1))

            for w in cands:
                in_k = w in R_k
                in_km1 = w in R_km1
                up, down = up_down(u, w, k, R_km1, N)

                if in_k and in_km1:  # Case 1 (Lemma 3.9)
                    counts["case1"] += 1
                    if up or down:
                        viol["case1_updown_nonempty"].append((tuple(u), k, tuple(w), len(up), len(down)))
                    # structural: u(k)=w(k), P_k = P_{k-1} u {u(k)}, n_k=n_{k-1}, q_k=q_{k-1}
                    if val(u, k) != val(w, k):
                        viol["case1_uk_ne_wk"].append((tuple(u), k, tuple(w)))
                    if P_set(u, w, k) != (P_set(u, w, k - 1) | {val(u, k)}):
                        viol["case1_Pk"].append((tuple(u), k, tuple(w)))
                    if n_val(u, w, k) != n_val(u, w, k - 1):
                        viol["case1_nk"].append((tuple(u), k, tuple(w)))
                    if R_k[w] != R_km1[w]:
                        viol["case1_qk"].append((tuple(u), k, tuple(w)))

                elif in_k and not in_km1:  # Case 2 (Lemma 3.8)
                    counts["case2"] += 1
                    if len(up) != 1 or len(down) != 0:
                        viol["case2_updown_shape"].append((tuple(u), k, tuple(w), len(up), len(down)))
                    if len(up) == 1:
                        b, t, v = up[0]
                        if P_set(u, w, k) != P_set(u, v, k - 1):
                            viol["case2_P"].append((tuple(u), k, tuple(w)))
                        if n_val(u, w, k) != n_val(u, v, k - 1) + 1:
                            viol["case2_n"].append((tuple(u), k, tuple(w)))
                        expected = R_km1[v] if t == "raise" else R_km1[v] * q_ab(k, b)
                        if R_k[w] != expected:
                            viol["case2_q"].append((tuple(u), k, tuple(w), t, str(R_k[w]), str(expected)))

                elif (not in_k) and in_km1:  # Case 3 (Lemma 3.10)
                    counts["case3"] += 1
                    if len(down) != 1 or len(up) != 0:
                        viol["case3_updown_shape"].append((tuple(u), k, tuple(w), len(up), len(down)))
                    if len(down) == 1:
                        a, t, v = down[0]
                        # paper claim: P_{k-1}(u,v) = P_{k-1}(u,w) u {u(k)}
                        if P_set(u, v, k - 1) != (P_set(u, w, k - 1) | {val(u, k)}):
                            viol["case3_P"].append((tuple(u), k, tuple(w), sorted(P_set(u, v, k - 1)), sorted(P_set(u, w, k - 1) | {val(u, k)})))

                else:  # Case 4 (Lemma 3.11)
                    counts["case4"] += 1
                    if len(up) != len(down) or len(up) > 1:
                        viol["case4_updown_shape"].append((tuple(u), k, tuple(w), len(up), len(down)))
                    if len(up) == 1 and len(down) == 1:
                        b, tb, vb = up[0]
                        a, ta, va = down[0]
                        if P_set(u, vb, k - 1) != P_set(u, va, k - 1):
                            viol["case4_P"].append((tuple(u), k, tuple(w)))
                        if n_val(u, vb, k - 1) != n_val(u, va, k - 1):
                            viol["case4_n"].append((tuple(u), k, tuple(w)))
                        wb = R_km1[vb] * (q_ab(k, b) if tb == "drop" else S.One)
                        wa = R_km1[va] * (q_ab(a, k) if ta == "drop" else S.One)
                        if wb != wa:
                            viol["case4_q"].append((tuple(u), k, tuple(w), str(wb), str(wa)))

    print(f"n={n}: case counts:", dict(counts))
    if not viol:
        print("  ALL per-case structural claims hold. (Lemmas 3.8-3.11 as used are correct.)")
    else:
        for kk, lst in viol.items():
            print(f"  VIOLATION [{kk}]: {len(lst)} instance(s); e.g. {lst[:3]}")
    return not viol


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    ok = run(n)
    sys.exit(0 if ok else 1)
