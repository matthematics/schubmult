"""This will be important!!!"""
from schubmult import *
from schubmult.symbolic.common_polys import grothendieck_poly
from schubmult.abc import x, y
from schubmult.rings.schubert.grothendieck_ring import GrothendieckRing

#probably need to go both ways
def hecke_pairs(w: Permutation):
    """Finds all u, v such that u @ v == w"""
    stack = [(w, Permutation([]))]
    ret = set()
    while len(stack) > 0:
        u, v = stack.pop()
        ret.add((u, v))

        iv = ~v
        
        for desc in u.descents():
            u2 = u.swap(desc, desc + 1)
            
            if iv[desc] < iv[desc + 1]:
                iv2 = iv.swap(desc, desc + 1)
                v2 = ~iv2
                if (u, v2) not in ret and (u@v2) == w:
                    stack.append((u, v2))
                if (u2, v2) not in ret and (u2@v2) == w:
                    stack.append((u2, v2))
            else:
                if (u2, v) not in ret and (u2@v) == w:
                    stack.append((u2, v))
                
    return ret

if __name__ == "__main__":
    import sys

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    Gy = GrothendieckRing(y)
    for perm in perms:
        if perm.inv == 0:
            continue
        perm_pairs = hecke_pairs(perm)
        
        # cauchy_groth = sum([Gx._beta**(p1.inv + p2.inv - perm.inv) * Gx(p2).expand() * Gy(~p1).expand() for p1, p2 in perm_pairs]).expand()
        # dgroth = grothendieck_poly(perm, x, y, Gx._beta).expand()
        # assert (cauchy_groth-dgroth).expand() == 0, f"Grothendieck Cauchy identity failed for {perm}: \n{cauchy_groth=}\n{dgroth=}\n{(cauchy_groth-dgroth).expand()}"
        # print("Hat buckets")

        cauchy_groth_tensor = sum([Gx._beta**(p1.inv + p2.inv - perm.inv) * Gx(p2) @ Gy(~p1) for p1, p2 in perm_pairs])

        if all(v == 1 for v in perm.trimcode):
            print(cauchy_groth_tensor)
            