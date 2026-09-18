from schubmult import *


if __name__ == "__main__":
    import sys
    import itertools
    from schubmult.combinatorics.crystal_graph import CrystalGraphTensor
    from schubmult.utils.tuple_utils import pad_tuple

    g = GrassTensorAlgebra()
    r = RCGraphRing()
    
    n = int(sys.argv[1])
    extra = 0
    perms = [perm for perm in Permutation.all_permutations(n + extra) if len(perm.trimcode) <= n]
    
    for perm1, perm2 in itertools.product(perms, repeat=2):
        if perm1.inv == 0 or perm2.inv == 0:
            continue
    
        result = r.zero

        prd = Sx(perm1) * Sx(perm2)
        for rc1, cem_dict in RCGraph.full_CEM(perm1, n).items():
            term_to_add = r.zero
            for rc2, cem_dict2 in RCGraph.full_CEM(perm2, n).items():
                
                g1 = g.from_dict(cem_dict)
                g2 = g.from_dict(cem_dict2)
                #graph_base = (g1 * g2).to_rc_graph_ring_element().resize(n)
                
                for tensor1, coeff1 in g1.items():
                    for tensor2, coeff2 in g2.items():
                        rc_base = next(iter((g(tensor1)* g(tensor2)).to_rc_graph_ring_element().resize(n)))
                        if rc_base.is_principal:
                            #assert coeff1 * coeff2 >= 0
                            result += coeff1 * coeff2 * r(rc_base)
        assert all(v >= 0 for v in result.values()), f"Failure for {perm1}, {perm2} with term {term_to_add}"
        
        
        prd2 = Sx.zero

        for rc, coeff in result.items():
            assert coeff == prd.get(rc.perm, 0), f"Failure for {perm1}, {perm2} at {rc.perm} with coeff {coeff} vs {prd.get(rc.perm, 0)}"
            prd2 += coeff * Sx(rc.perm)
        
        assert prd == prd2, f"Failure for {perm1}, {perm2} with result {result}, expected {prd}, got {prd2}"
        print("Success for", perm1, perm2)

        
        
