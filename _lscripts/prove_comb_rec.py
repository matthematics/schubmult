"""Attack COMB-REC directly by the position-k structural axis.

COMB-REC (to prove) for all u,w and 1<=p<=k:
  c_k(p;u,w) = (y_{w(k)}-z_j) c_{k-1}(p-1;u,w)                                  [DIAG]
      + sum_{pp!=k} sign(pp-k)[q_{k,pp} if drop] c_{k-1}(p-1; u, w t_{k,pp})    [MONK]
      + c_{k-1}(p;u,w)                                                          [REC]
      + q_{k-1} c_{k-2}(p-2;u,w)                                                [QUANT]
where j=k-p+1 and c_m(a;u,w)=q_m(u,w) E_{a-n_m,m-n_m}(y_{P_m};z).

Hand result:
  Case A (u(k)=w(k)):  DIAG+REC == c_k(p;u,w) EXACTLY (factorial-elem single-var recursion,
      adding variable y_{w(k)}).  => residual R_A := MONK + QUANT must be 0.
  Case B (u(k)!=w(k)): c_k(p;u,w) = (q_k/q_{k-1}-shifted) c_{k-1}(p-1;u,w)-type term; DIAG,REC
      restructure differently.  => residual R_B := c_k - (DIAG+MONK+REC+QUANT) must be 0
      but we isolate which pieces carry it.

This script:
  (1) confirms Case A: DIAG+REC == c_k  and  R_A == 0, separately, for all fixed-k (u,w).
  (2) confirms Case B: isolates the identity in the non-fixed case.
It PRINTS per-case pass/fail so we know exactly which sub-identity to prove.

Run: conda activate schubmult_312 && python _lscripts/prove_comb_rec.py 4
"""

import sys
from collections import defaultdict

from schubmult import *  # noqa: F401,F403
from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import S, expand, expand_func, prod
from schubmult.symbolic.poly.variables import GeneratingSet
from schubmult.symbolic.symmetric_polynomials import FactorialElemSym

y = GeneratingSet("y")
z = GeneratingSet("z")
q_gs = GeneratingSet("q")


def q_ab(a, b):
    return prod([q_gs[s] for s in range(a, b)])


def val(perm, pos):
    return perm[pos - 1]


def enumerate_pieri(u, k, N):
    """{w: q_k(u,w)} for all w with u ->_k w (Definition 1.1)."""
    u = Permutation(u)
    results = {u: S.One}
    stack = [(u, frozenset(), N + 1, S.One, u.inv)]
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
                results.setdefault(nperm, nqw)
                stack.append((nperm, used_a | {a}, b, nqw, nperm.inv))
    return results


def c_term(u, w, m, a, Rm):
    """c_m(a;u,w) = q_m E_{a-n_m, m-n_m}(y_{P_m}; z), 0 unless u->_m w and n_m<=a<=m."""
    if m < 0:
        return S.Zero
    w = Permutation(w)
    if w not in Rm:
        return S.Zero
    qm = Rm[w]
    Pm = sorted(val(u, i) for i in range(1, m + 1) if val(u, i) == val(w, i))
    nm = m - len(Pm)
    deg = a - nm
    nv = m - nm
    if deg < 0 or deg > nv:
        return S.Zero
    if nv == 0:
        return qm if deg == 0 else S.Zero
    yvars = [y[v] for v in Pm]
    ncoeff = nv + 1 - deg
    zvars = [z[i] for i in range(1, ncoeff + 1)]
    return qm * expand_func(FactorialElemSym(deg, nv, yvars, zvars))


def diag_term(u, p, k, w, Rkm1):
    j = k - p + 1
    return (y[val(w, k)] - z[j]) * c_term(u, w, k - 1, p - 1, Rkm1)


def monk_term(u, p, k, w, N, Rkm1):
    total = S.Zero
    for pp in range(1, N + 1):
        if pp == k:
            continue
        v = w.swap(k - 1, pp - 1)
        dl = w.inv - v.inv
        sign = 1 if pp > k else -1
        lo, hi = min(k, pp), max(k, pp)
        if dl == 1:
            total += sign * c_term(u, v, k - 1, p - 1, Rkm1)
        elif dl == 1 - 2 * (hi - lo):
            total += sign * q_ab(lo, hi) * c_term(u, v, k - 1, p - 1, Rkm1)
    return total


def rec_term(u, p, k, w, Rkm1):
    return c_term(u, w, k - 1, p, Rkm1)


def quant_term(u, p, k, w, Rkm2):
    return q_gs[k - 1] * c_term(u, w, k - 2, p - 2, Rkm2)


def run(n):
    N = 2 * n
    perms = list(Permutation.all_permutations(n))
    stats = defaultdict(lambda: {"checks": 0, "fails": 0, "examples": []})
    for u in perms:
        u = Permutation(u)
        for k in range(2, n):
            Rk = enumerate_pieri(u, k, N)
            Rkm1 = enumerate_pieri(u, k - 1, N)
            Rkm2 = enumerate_pieri(u, k - 2, N) if k - 2 >= 0 else {u: S.One}
            cands = set(Rk) | set(Rkm1) | set(Rkm2)
            for v in list(Rkm1):
                for pp in range(1, N + 1):
                    if pp != k:
                        cands.add(v.swap(k - 1, pp - 1))
            for p in range(1, k + 1):
                for w in cands:
                    w = Permutation(w)
                    fixed_k = val(u, k) == val(w, k)
                    ck = c_term(u, w, k, p, Rk)
                    DIAG = diag_term(u, p, k, w, Rkm1)
                    MONK = monk_term(u, p, k, w, N, Rkm1)
                    REC = rec_term(u, p, k, w, Rkm1)
                    QUANT = quant_term(u, p, k, w, Rkm2)
                    case = "A_fixed" if fixed_k else "B_moved"
                    # Full COMB-REC (sanity)
                    s = stats[(case, "FULL")]
                    s["checks"] += 1
                    if expand(ck - (DIAG + MONK + REC + QUANT)) != 0:
                        s["fails"] += 1
                        if len(s["examples"]) < 5:
                            s["examples"].append((tuple(u), p, k, tuple(w)))
                    if fixed_k:
                        # Case A sub-identities
                        s1 = stats[("A_fixed", "DIAG+REC==ck")]
                        s1["checks"] += 1
                        if expand(ck - (DIAG + REC)) != 0:
                            s1["fails"] += 1
                            if len(s1["examples"]) < 5:
                                s1["examples"].append((tuple(u), p, k, tuple(w)))
                        s2 = stats[("A_fixed", "R_A:MONK+QUANT==0")]
                        s2["checks"] += 1
                        if expand(MONK + QUANT) != 0:
                            s2["fails"] += 1
                            if len(s2["examples"]) < 5:
                                s2["examples"].append((tuple(u), p, k, tuple(w)))
                    else:
                        in_km1 = Permutation(w) in Rkm1
                        if not in_km1:
                            # Case B2: u ->_k w but NOT u ->_{k-1} w. DIAG=REC=0.
                            # Claim: c_k == MONK + QUANT.
                            s3 = stats[("B_moved", "B2(no km1): ck==MONK+QUANT")]
                            s3["checks"] += 1
                            if expand(ck - (MONK + QUANT)) != 0:
                                s3["fails"] += 1
                                if len(s3["examples"]) < 5:
                                    s3["examples"].append((tuple(u), p, k, tuple(w)))
                        else:
                            # Case B1: u ->_k w AND u ->_{k-1} w.
                            # Full identity; isolate MONK+QUANT residual.
                            s4 = stats[("B_moved", "B1(yes km1): ck-DIAG-REC==MONK+QUANT")]
                            s4["checks"] += 1
                            if expand((ck - DIAG - REC) - (MONK + QUANT)) != 0:
                                s4["fails"] += 1
                                if len(s4["examples"]) < 5:
                                    s4["examples"].append((tuple(u), p, k, tuple(w)))
    print(f"n={n}:")
    for (case, name), s in sorted(stats.items()):
        flag = "OK" if s["fails"] == 0 else f"FAIL({s['fails']})"
        print(f"  [{case}] {name}: {s['checks']} checks -> {flag}")
        for ex in s["examples"]:
            print(f"      ex: u={ex[0]} p={ex[1]} k={ex[2]} w={ex[3]}")
    return all(s["fails"] == 0 for s in stats.values())


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    ok = run(n)
    print("RESULT:", "all sub-identities hold" if ok else "some sub-identity fails")
    sys.exit(0 if ok else 1)
