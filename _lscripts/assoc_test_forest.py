from schubmult import *
from schubmult.rings.combinatorial.forest_rc_ring import ForestRCGraphRing, DualForestRCGraphRing

def the_forests(perm, length):
    return {rc for rc in RCGraph.all_rc_graphs(perm, length) if rc.forest_weight == perm.pad_code(length)}

if __name__ == "__main__":
    import sys
    import itertools
    n = int(sys.argv[1])
    r = ForestRCGraphRing()
    dr = DualForestRCGraphRing()
    perms = Permutation.all_permutations(n)
    failures = {}
    length = n - 1
    for perm1 in perms:
        for rc1 in the_forests(perm1, length):
            for perm2 in perms:
                for rc2 in the_forests(perm2, length): 
                    for perm3 in perms:
                        for rc3 in the_forests(perm3, length):
                            print(f"Checking {repr(rc1)}, {repr(rc2)}, {repr(rc3)}")
                            snap_size1 = max(len(rc1.perm.trimcode), len(rc2.perm.trimcode))
                            snap_size2 = max(snap_size1, len(rc3.perm.trimcode))
                            snap_size3 = max(len(rc2.perm.trimcode), len(rc3.perm.trimcode))
                            
                            prd = (r(rc1) *  (r(rc2) *  r(rc3)))
                            prd2 = ((r(rc1) *  r(rc2)) *  r(rc3))
                            if not prd.almosteq(prd2):
                                if not prd.almosteq(0) and not prd2.almosteq(0):
                                    print(f"Product mismatch for {rc1}, {rc2}, {rc3}: expected {prd}, got {prd2}")
                                    failures[(rc1, rc2, rc3)] = (prd, prd2)
                            else:
                                print(f"Products match")
    # for perm1 in perms:
    #     for rc1 in the_forests(perm1, length):
    #         for perm2 in perms:
    #             for rc2 in the_forests(perm2, length): 
    #         #         for perm3 in perms:
    #         #             for rc3 in the_forests(perm3, length):
    #         #                 print(f"Checking {repr(rc1)}, {repr(rc2)}, {repr(rc3)}")
                            
    #         #                 prd = (dr(rc1) *  (dr(rc2) *  dr(rc3)))
    #         #                 prd2 = ((dr(rc1) *  dr(rc2)) *  dr(rc3))
    #         #                 if not prd.almosteq(prd2):
    #         #                     if not prd.almosteq(0) and not prd2.almosteq(0):
    #         #                         print(f"Product mismatch for {rc1}, {rc2}, {rc3}: expected {prd-prd2}")
    #         #                         failures[(rc1, rc2, rc3)] = (prd, prd2)
    #         #                 else:
    #         #                     print(f"Products match")
    #                 print(f"Checking coproduct for {repr(rc1)}, {repr(rc2)}")
    #                 prd = r(rc1).coproduct() * r(rc2).coproduct()
    #                 prd2 = (r(rc1) * r(rc2)).coproduct()
    #                 if not prd.almosteq(prd2):
    #                     print(f"Product mismatch for {rc1}, {rc2}: expected: {prd} \ngot: {prd2}")
    #                     failures[(rc1, rc2)] = (prd, prd2)
    if failures:
        print(f"Failures: {len(failures)}")
        for (rc1, rc2), (prd, prd2) in failures.items():
            print(f"Failure: {repr(rc1)}, {repr(rc2)}")
            print(f"Expected: {prd}, Got: {prd2}")
    else:
        print("All tests passed!")
                        
