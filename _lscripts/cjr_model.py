"""FTFSDL / canonical join complex model for counting Lambda.

Claims to test on weak order (semidistributive):
 (A) x <= w  iff  CJR(x) subset J(w),  where J(w) = {join-irreducibles <= w}.
     Hence Lambda(w) = # faces of the canonical join complex supported on J(w).
 (B) the canonical join complex is FLAG (Barnard): a set of join-irreducibles is a
     face iff pairwise compatible.  So Lambda(w) = # cliques in a graph.
 (C) split join-irreducibles into U (uncontracted by sylvester) and C (contracted).
     Every clique K splits as (K∩U, K∩C), injectively, so
        Lambda(w) <= #cliques(G_w|U) * #cliques(G_w|C).
     Is #cliques(G_w|U) = #dominants <= w, and #cliques(G_w|C) related to class size?
"""

import sys
from itertools import combinations

from schubmult import Permutation, uncode


def main(n):
    P = list(Permutation.all_permutations(n))
    inv = {w: frozenset(w.inversion_set) for w in P}

    def leq(a, b):
        return inv[a] <= inv[b]

    # covers
    lower = {w: [] for w in P}
    for x in P:
        for y in P:
            if leq(y, x) and len(inv[y]) == len(inv[x]) - 1:
                lower[x].append(y)

    JI = [x for x in P if len(lower[x]) == 1]
    print(f"S_{n}: {len(P)} elements, {len(JI)} join-irreducibles")

    # canonical join representation
    def cjr(x):
        out = []
        for y in lower[x]:
            cands = [z for z in P if leq(z, x) and not leq(z, y)]
            mn = min(cands, key=lambda z: len(inv[z]))
            assert all(leq(mn, z) for z in cands), f"no unique minimum for {x},{y}"
            out.append(mn)
        return frozenset(out)

    CJR = {x: cjr(x) for x in P}
    assert all(all(j in JI for j in CJR[x]) for x in P), "CJR not made of join-irreducibles"

    # (A)
    badA = 0
    for w in P:
        Jw = {j for j in JI if leq(j, w)}
        for x in P:
            if leq(x, w) != (CJR[x] <= Jw):
                badA += 1
    print(f"  (A) x<=w  iff  CJR(x) subset J(w) : {badA} failures")

    # bijection onto faces
    faces = {CJR[x] for x in P}
    print(f"  CJR injective: {len(faces) == len(P)}")

    # (B) flagness
    compat = {frozenset(f) for f in faces if len(f) == 2}
    badB = 0
    for x in P:
        F = list(CJR[x])
        if not all(frozenset({a, b}) in compat for a, b in combinations(F, 2)):
            badB += 1
    pairs_ok = 0
    for a, b in combinations(JI, 2):
        if frozenset({a, b}) in compat and frozenset({a, b}) not in faces:
            pairs_ok += 1
    # flag: every pairwise-compatible set is a face
    badflag = 0
    allsets = 0
    for k in range(1, min(len(JI), 5) + 1):
        for S in combinations(JI, k):
            if all(frozenset({a, b}) in compat for a, b in combinations(S, 2)):
                allsets += 1
                if frozenset(S) not in faces:
                    badflag += 1
    print(f"  (B) every pairwise-compatible set (size<=5) is a face: {badflag} failures / {allsets}")

    # (C) split by sylvester congruence
    def mdom(v):
        return Permutation(list(~uncode((~v).theta())))

    U = [j for j in JI if mdom(j) == j]   # uncontracted: j is its own class top
    Cc = [j for j in JI if mdom(j) != j]
    print(f"  join-irreducibles: {len(U)} with mdom(j)==j, {len(Cc)} others")

    def ncliques(verts):
        vs = list(verts)
        cnt = 0
        for k in range(len(vs) + 1):
            for S in combinations(vs, k):
                if all(frozenset({a, b}) in compat for a, b in combinations(S, 2)):
                    cnt += 1
        return cnt

    badC = 0
    worst = None
    doms = [w for w in P if mdom(w) == w]
    classes = {}
    for v in P:
        classes.setdefault(mdom(v), []).append(v)
    for w in doms:
        Jw = [j for j in JI if leq(j, w)]
        lam = sum(1 for x in P if leq(x, w))
        cu = ncliques([j for j in Jw if j in U])
        cc = ncliques([j for j in Jw if j in Cc])
        ndom = sum(1 for d in doms if leq(d, w))
        m = len(classes[w])
        if lam > cu * cc:
            badC += 1
        r = lam / (cu * cc)
        if worst is None or r < worst[0]:
            worst = (r, w, lam, cu, cc, ndom, m)
    print(f"  (C) Lambda(w) <= cliques(U)*cliques(C): {badC} failures over {len(doms)} dominants")
    r, w, lam, cu, cc, ndom, m = worst
    print(f"      loosest: w={w} Lambda={lam} cliquesU={cu} cliquesC={cc} #dom<=w={ndom} |class|={m}")

    # is cliques(U) exactly the number of dominants below w?
    badU = 0
    for w in doms:
        Jw = [j for j in JI if leq(j, w)]
        cu = ncliques([j for j in Jw if j in U])
        ndom = sum(1 for d in doms if leq(d, w))
        if cu != ndom:
            badU += 1
    print(f"  cliques(U-part) == #dominants below w : {badU} failures")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 4)
