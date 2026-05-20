from schubmult import *
from schubmult.rings.combinatorial.slide_rc_ring import SlideRCGraphRing

def the_slides(perm, length):
    return {rc for rc in RCGraph.all_rc_graphs(perm, length) if rc.snap_qy().length_vector == perm.pad_code(length)}

if __name__ == "__main__":
    import sys
    import itertools
    n = int(sys.argv[1])
    r = SlideRCGraphRing()
    perms = [perm for perm in Permutation.all_permutations(n) if perm.inv > 0]
    failures = {}
    length = n - 1
    for perm1 in perms:
        for rc1 in the_slides(perm1, length):
            for perm2 in perms:
                for rc2 in the_slides(perm2, length): 
                    for perm3 in perms:
                        for rc3 in the_slides(perm3, length):
                            length0 = max(len(perm1.trimcode), len(perm2.trimcode), len(perm3.trimcode))
                            #length0 = length
                            rc1_sized = rc1.resize(length0)
                            rc2_sized = rc2.resize(length0)
                            rc3_sized = rc3.resize(length0)
                            print(f"Checking {repr(rc1_sized)}, {repr(rc2_sized)}, {repr(rc3_sized)}")
                            
                            prd = (r(rc1_sized) * r(rc2_sized)) * r(rc3_sized)
                            prd2 = (r(rc1_sized) * (r(rc2_sized) * r(rc3_sized)))
                            if prd == 0:
                                prd = r.zero
                            if not prd.almosteq(prd2):
                                print(f"Product mismatch for {rc1}, {rc2}, {rc3}: expected {prd}, got {prd2}")
                                failures[(rc1, rc2, rc3)] = (prd, prd2)
                                raise ValueError(f"Product mismatch for {rc1}, {rc2}, {rc3}: expected {prd}, got {prd2}")
                            else:
                                print(f"Products match")
    if failures:
        print(f"Failures: {len(failures)}")
        for (rc1, rc2, rc3), (prd, prd2) in failures.items():
            print(f"Failure: {repr(rc1)}, {repr(rc2)}, {repr(rc3)}")
            print(f"Expected: {prd}, Got: {prd2}")
    else:
        print("All tests passed!")
                        
