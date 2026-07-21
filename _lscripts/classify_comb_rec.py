"""Decisive classification of COMB-REC into the corrected, minimal proof skeleton.

For every (u,p,k,w), classify by the triple
    (fixed = [u(k)==w(k)],  ink = [u->_k w],  inkm1 = [u->_{k-1} w])
and verify the reduced identity that holds in each occurring class.  Also verify the
sharp structural facts that make R_A and the B1 residual PROVABLE:

  A (fixed):        DIAG+REC == c_k   AND   R_A: MONK+QUANT == 0
                    sharp: MONK has <=1 nonzero term, a down-neighbor, == -QUANT.
  B2 (moved, ink, not inkm1):  c_k == MONK+QUANT     (= paper Lemma 3.8, proven)
  B0 (moved, not ink):         c_k = 0;  verify DIAG+MONK+REC+QUANT == 0 too.
  B1 (moved, ink, inkm1):      does this class even occur?  If empty, B1 is vacuous.

Prints the population of every triple and per-class pass/fail.

Run: conda activate schubmult_312 && python _lscripts/classify_comb_rec.py 5
"""

import sys
from collections import Counter, defaultdict

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


def monk_terms_list(u, p, k, w, N, Rkm1):
    """List of (pp, dir, type, term) for nonzero Monk off-diagonal contributions."""
    out = []
    for pp in range(1, N + 1):
        if pp == k:
            continue
        vv = w.swap(k - 1, pp - 1)
        dl = w.inv - vv.inv
        sign = 1 if pp > k else -1
        lo, hi = min(k, pp), max(k, pp)
        if dl == 1:
            term = sign * c_term(u, vv, k - 1, p - 1, Rkm1)
            ttype = "RAISE"
        elif dl == 1 - 2 * (hi - lo):
            term = sign * q_ab(lo, hi) * c_term(u, vv, k - 1, p - 1, Rkm1)
            ttype = "DROP"
        else:
            continue
        if expand(term) != 0:
            out.append((pp, "up" if pp > k else "dn", ttype, term))
    return out


def run(n):
    N = 2 * n
    perms = list(Permutation.all_permutations(n))
    triple_pop = Counter()
    class_stats = defaultdict(lambda: {"checks": 0, "fails": 0, "ex": []})
    ra_monk_hist = Counter()
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
                    fixed = val(u, k) == val(w, k)
                    ink = w in Rk
                    inkm1 = w in Rkm1
                    triple_pop[(fixed, ink, inkm1)] += 1
                    j = k - p + 1
                    DIAG = (y[val(w, k)] - z[j]) * c_term(u, w, k - 1, p - 1, Rkm1)
                    mlist = monk_terms_list(u, p, k, w, N, Rkm1)
                    MONK = sum((t[3] for t in mlist), S.Zero)
                    REC = c_term(u, w, k - 1, p, Rkm1)
                    QUANT = q_gs[k - 1] * c_term(u, w, k - 2, p - 2, Rkm2)
                    ck = c_term(u, w, k, p, Rk)

                    def check(label, expr):
                        st = class_stats[label]
                        st["checks"] += 1
                        if expand(expr) != 0:
                            st["fails"] += 1
                            if len(st["ex"]) < 4:
                                st["ex"].append((tuple(u), p, k, tuple(w)))

                    if fixed:
                        check("A: DIAG+REC==ck", ck - (DIAG + REC))
                        check("A: R_A MONK+QUANT==0", MONK + QUANT)
                        # sharp: down-only, <=1 term
                        dn = [t for t in mlist if t[1] == "dn"]
                        up = [t for t in mlist if t[1] == "up"]
                        ra_monk_hist[(len(dn), len(up))] += 1
                    else:
                        if ink and not inkm1:
                            check("B2: ck==MONK+QUANT", ck - (MONK + QUANT))
                        elif not ink:
                            check("B0(ck=0): full==0", DIAG + MONK + REC + QUANT)
                        else:  # ink and inkm1
                            check("B1: ck-DIAG-REC==MONK+QUANT", (ck - DIAG - REC) - (MONK + QUANT))
    print(f"n={n}: triple population (fixed,ink,inkm1):")
    for key, cnt in sorted(triple_pop.items()):
        print(f"    {key}: {cnt}")
    print("  class identities:")
    for label, st in sorted(class_stats.items()):
        flag = "OK" if st["fails"] == 0 else f"FAIL({st['fails']})"
        print(f"    {label}: {st['checks']} -> {flag}")
        for e in st["ex"]:
            print(f"        ex {e}")
    print(f"  R_A sharp (dn,up) MONK-term-count histogram: {dict(sorted(ra_monk_hist.items()))}")


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    run(n)
