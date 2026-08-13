"""Lambda(w) = #cliques in the arc-compatibility graph induced on J(w).
Uses recursive clique counting (cheap: total cliques = |L|).
Also retests the U/C split with the CORRECT contracted definition:
   j is contracted by a congruence iff j == j_* modulo the congruence.
"""
import sys
from itertools import combinations
from schubmult import Permutation, uncode


def main(n):
    P = list(Permutation.all_permutations(n))
    inv = {w: frozenset(w.inversion_set) for w in P}
    def leq(a, b): return inv[a] <= inv[b]

    lower = {w: [] for w in P}
    for x in P:
        lx = len(inv[x])
        for y in P:
            if len(inv[y]) == lx - 1 and leq(y, x):
                lower[x].append(y)
    JI = [x for x in P if len(lower[x]) == 1]
    jstar = {j: lower[j][0] for j in JI}

    def cjr(x):
        out = []
        for y in lower[x]:
            cands = [z for z in P if leq(z, x) and not leq(z, y)]
            mn = min(cands, key=lambda z: len(inv[z]))
            assert all(leq(mn, z) for z in cands)
            out.append(mn)
        return frozenset(out)
    CJR = {x: cjr(x) for x in P}

    badA = 0
    for w in P:
        Jw = {j for j in JI if leq(j, w)}
        for x in P:
            if leq(x, w) != (CJR[x] <= Jw):
                badA += 1
    faces = {CJR[x] for x in P}
    compat = {f for f in faces if len(f) == 2}
    print(f"S_{n}: {len(P)} elements, {len(JI)} join-irreducibles")
    print(f"  (A) x<=w iff CJR(x) subset J(w): {badA} failures")
    print(f"  CJR injective onto faces: {len(faces) == len(P)}")

    badflag = 0
    for x in P:
        F = list(CJR[x])
        if not all(frozenset({a, b}) in compat for a, b in combinations(F, 2)):
            badflag += 1
    print(f"  (B) faces are pairwise-compatible: {badflag} failures")

    adj = {j: {k for k in JI if k != j and frozenset({j, k}) in compat} for j in JI}

    def count_cliques(verts):
        vs = sorted(verts, key=lambda z: (len(inv[z]), tuple(z)))
        idx = {v: i for i, v in enumerate(vs)}
        total = 0
        def rec(cur_i, cand):
            nonlocal total
            total += 1
            for v in list(cand):
                if idx[v] < cur_i:
                    continue
                rec(idx[v] + 1, cand & adj[v])
        rec(0, set(vs))
        return total

    badC = 0
    for w in P:
        Jw = [j for j in JI if leq(j, w)]
        lam = sum(1 for x in P if leq(x, w))
        if count_cliques(Jw) != lam:
            badC += 1
    print(f"  Lambda(w) == #cliques on J(w): {badC} failures over {len(P)} elements")

    def mdom(v): return Permutation(list(~uncode((~v).theta())))
    contracted = [j for j in JI if mdom(j) == mdom(jstar[j])]
    unc = [j for j in JI if mdom(j) != mdom(jstar[j])]
    doms = [w for w in P if mdom(w) == w]
    print(f"  join-irreducibles: {len(unc)} uncontracted, {len(contracted)} contracted"
          f"  (Tamari has {len(doms)} elements, expect n(n-1)/2={n*(n-1)//2} uncontracted)")

    badU = 0
    worst = None
    classes = {}
    for v in P:
        classes.setdefault(mdom(v), []).append(v)
    for w in doms:
        Jw = [j for j in JI if leq(j, w)]
        lam = sum(1 for x in P if leq(x, w))
        cu = count_cliques([j for j in Jw if j in unc])
        cc = count_cliques([j for j in Jw if j in contracted])
        nd = sum(1 for d in doms if leq(d, w))
        if cu != nd:
            badU += 1
        r = lam / (cu * cc) if cu * cc else 0
        if worst is None or r > worst[0]:
            worst = (r, w, lam, cu, cc, nd, len(classes[w]))
    print(f"  #cliques(uncontracted part) == #dominants below w: {badU} failures / {len(doms)}")
    r, w, lam, cu, cc, nd, m = worst
    print(f"     tightest split: w={w} Lambda={lam} cu={cu} cc={cc} #dom<=w={nd} |class|={m} ratio={r:.3f}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
