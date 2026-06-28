from schubmult import *
from schubmult.symbolic.poly.variables import MaskedGeneratingSet
from schubmult.abc import x
from schubmult.symbolic import Symbol

if __name__ == "__main__":
    import sys

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n + 1)
    b = Symbol("b")
    #for perm in perms:
    results = {}
    mu2 = ~uncode(list(reversed(range(1, (n // 2) + 1))))
    mu1 = ~uncode(list(reversed(range((n // 2) + 1, n + 1))))
    gset2 = MaskedGeneratingSet(x, list(range(1, (n-(n//2)) + 1)) + list(range(20,50)))
    gset1 = gset2.complement()
    print(gset1[:len(gset1)], gset2[:10])
    Gx1 = GrothendieckRing(gset1)
    Gx2 = GrothendieckRing(gset2)
    for u1 in perms:
        for v1 in perms:
            w1 = (~u1) @ v1
            if w1 != mu1:
                continue
            for u2 in perms:
                for v2 in perms:
                    w2 = (~u2) @ v2
                    if w2 != mu2:
                        continue
                    # if u1.inv + v1.inv - w1.inv + u2.inv + v2.inv - w2.inv != 0:
                    #     continue
                    the_prod = b ** (u1.inv + v1.inv - w1.inv + u2.inv + v2.inv - w2.inv) * Gx(v1) * Gx(v2)
                    for w, coeff in the_prod.items():
                        results[w] = results.get(w, 0) + coeff * (Gx1(u1)@Gx2(u2))
    for w, coeff in results.items():
        pol1 = Gx(w*Permutation.w0(n+1)).as_polynomial()
        pol2 = 0
        for (u1, u2), coeff2 in coeff.items():
            if coeff2 != 0:
                pol2 += coeff2 * Gx1(u1).as_polynomial() * Gx2(u2).as_polynomial()
        assert (pol1 - pol2.subs(b, 0)).expand() == 0, f"Mismatch for {w}: {(pol1-pol2.subs(b, 0)).expand()}\n{pol1}\n{pol2.subs(b,0)}"
        print(f"Verified {w}")