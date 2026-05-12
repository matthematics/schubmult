from schubmult import *
from schubmult.rings.free_algebra import *
from schubmult.rings.combinatorial.qy_rc_graph_ring import _canonical_rc, QYRCGraphRing
from schubmult.utils.perm_utils import weak_compositions, artin_sequences

r = RCGraphRing()
q = QYRCGraphRing()
FSlideDual = FreeAlgebra(FundamentalSlideBasis)
def rc_slide_elem(comp):
    return sum([r(rc) for rc in r.monomial(*comp).keys() if rc.is_quasi_yamanouchi])

def slide_coproduct_rc(comp):
    cprd = FSlideDual(*comp).coproduct()
    result = (r@r).zero
    for (comp1, comp2), coeff in cprd.items():
        result += coeff * rc_slide_elem(comp1) @ rc_slide_elem(comp2)
    return result

def verify_slide_coproduct(comp):
    cprd = slide_coproduct_rc(comp)
    result = (r@r).zero
    for (rc1, rc2), coeff in cprd.items():
        result += coeff * r(rc1) @ r(rc2)
    return result



def slide_hom(rc):
    if isinstance(rc, RCGraph):
        return q(_canonical_rc(rc)) @ FA(*rc.length_vector)
    result = 0
    for rc, coeff in rc.items():
        result += coeff * q(_canonical_rc(rc)) @ FA(*rc.length_vector)
    return result

def verify_slide_hom(n):
    import itertools
    
    perms = Permutation.all_permutations(n)
    print("Testing slide homomorphism for all pairs of RC graphs of permutations in S_{}...".format(n))
    for perm1, perm2 in itertools.product(perms, repeat=2):
        for rc1, rc2  in itertools.product(RCGraph.all_rc_graphs(perm1), RCGraph.all_rc_graphs(perm2)):
            test1 = slide_hom(rc1) * slide_hom(rc2)
            the_tester = slide_hom(r(rc1) * r(rc2))
            assert test1.almosteq(the_tester), f"Failed for {rc1}, {rc2} with test1={test1} and tester={the_tester}"

if __name__ == "__main__":
    import sys
    from sympy import pretty_print
    n = int(sys.argv[1])
    #max_degree = int(sys.argv[2])
    verify_slide_hom(n)
    
    
    