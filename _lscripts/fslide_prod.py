from schubmult import *
from schubmult.rings.polynomial_algebra import *
from schubmult.rings.free_algebra import *
from schubmult.rings.combinatorial.qy_rc_graph_ring import QYRCGraphRing, _canonical_rc
from sympy import pretty_print

r = RCGraphRing()
q = QYRCGraphRing()

def hom(rc_elem):
    result = 0
    for rc, coeff in rc_elem.items():
        result += coeff * (q(_canonical_rc(rc)) @ FA(*rc.length_vector))
    return result

def inv_hom(qyw_elem):
    result = 0
    for (rc, weight), coeff in qyw_elem.items():
        wrd = rc.perm_word
        comp_seq = []
        for i, val  in enumerate(weight, start=1):
            if val != 0:
                comp_seq.extend([i] * val)
        result += coeff * r(RCGraph.from_reduced_compatible(wrd, tuple(comp_seq)).resize(len(rc)))
    return result

if __name__ == "__main__":
    import sys
    import itertools

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)

    
    #FSlideDual = FreeAlgebra(FundamentalSlideBasis)

    for perm1, perm2 in itertools.product(perms, repeat=2):
        # comp1 = perm1.trimcode
        # comp2 = perm2.trimcode
        
        for length1, length2 in itertools.product(range(max(1,len(perm1.trimcode)), n), range(max(1,len(perm2.trimcode)), n)):
            for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1, length1), RCGraph.all_rc_graphs(perm2, length2)):
                test0 = r(rc1) * r(rc2)
                test = hom(test0)
                toast = hom(r(rc1)) * hom(r(rc2))
                assert test.almosteq(toast), f"Failed for {rc1=} {rc2=} {test=} {toast=}"
                assert inv_hom(test).almosteq(test0), f"Failed for {rc1=} {rc2=} {test=} {inv_hom(test)=}"