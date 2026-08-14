from schubmult import *
import sympy

r = RCGraphRing()

def fa_hom(rc):
    if isinstance(rc, RCGraph):
        return FA(*rc.length_vector)
    res = FA.zero
    for rc1, coeff in rc.items():
        res += coeff * fa_hom(rc1)
    return res

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    
    
    for perm1 in perms:
        for perm2 in perms:
            for k1 in range(max(1,perm1.max_descent), n):
                for rc1 in RCGraph.all_rc_graphs(perm1, k1):
                    for k2 in range(max(1,perm2.max_descent), n):
                        for rc2 in RCGraph.all_rc_graphs(perm2, k2):
                            prod = r.monomial(*rc1.length_vector) * r.monomial(*rc2.length_vector)
                            fa_prod = fa_hom(r.monomial(*rc1.length_vector)) * fa_hom(r.monomial(*rc2.length_vector))
                            test_fa_prod = fa_hom(prod)
                            assert test_fa_prod.almosteq(fa_prod), f"Discrepancy for {rc1} * {rc2}: {test_fa_prod} != {fa_prod}"
                

                