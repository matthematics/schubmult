from schubmult import *
from schubmult.rings.combinatorial.slide_rc_ring import SlideRCGraphRing
from schubmult.rings.combinatorial.schubert_rc_ring import SchubertRCGraphRing

#r = SlideRCGraphRing()
r = SchubertRCGraphRing()
br = BoundedRCFactorAlgebra()

def the_slides(perm, length):
    ret = {}
    for rc in RCGraph.all_rc_graphs(perm, length):
        #ret[rc.perm_word] = ret.get(rc.perm_word, 0) + r(rc)
        if rc.snap_qy().length_vector == perm.pad_code(length):
            ret[rc] = rc
    #ret[1] = r.from_dict(dict.fromkeys(RCGraph.all_rc_graphs(perm, length),1))
    return ret

if __name__ == "__main__":
    import sys
    import itertools
    n = int(sys.argv[1])
    #r = SlideRCGraphRing()
    perms = [perm for perm in Permutation.all_permutations(n) if perm.inv > 0]
    failures = {}
    length = n - 1
    for perm1 in perms:
        for word1, rc1 in the_slides(perm1, length).items():
            rc1_factor = br.from_rc_graph_ring_element(r(rc1), size=len(perm1))
            for perm2 in perms:
                for word2, rc2 in the_slides(perm2, length).items(): 
                    rc2_factor = br.from_rc_graph_ring_element(r(rc2), size=len(perm2))
                    for perm3 in perms:
                        for word3, rc3 in the_slides(perm3, length).items():
                            rc3_factor = br.from_rc_graph_ring_element(r(rc3), size=len(perm3))
                            prd = ((rc1_factor * rc2_factor) * rc3_factor).to_rc_graph_ring_element()
                            prd2 = (rc1_factor * (rc2_factor * rc3_factor)).to_rc_graph_ring_element()
                            if prd == 0:
                                prd = r.zero
                            if not prd.almosteq(prd2):
                                print(f"Product mismatch for {rc1}, {rc2}, {rc3}: expected {prd}, got {prd2}")
                                failures[(rc1, rc2, rc3)] = (prd, prd2)
                                raise ValueError(f"Product mismatch for {rc1}, {rc2}, {rc3}: expected {prd}, got {prd2}")
                            else:
                                print(f"Products match {perm1} {perm2} {perm3}")
    if failures:
        print(f"Failures: {len(failures)}")
        for (rc1, rc2, rc3), (prd, prd2) in failures.items():
            print(f"Failure: {repr(rc1)}, {repr(rc2)}, {repr(rc3)}")
            print(f"Expected: {prd}, Got: {prd2}")
    else:
        print("All tests passed!")
                        
