from schubmult import *
from schubmult.rings.polynomial_algebra import *
from schubmult.rings.free_algebra import *
from sympy import pretty_print
import itertools

def schub_coeff(perm, length_vec):
    return ASx(perm, len(length_vec)).change_basis(WordBasis).get(length_vec, 0)

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)
    
    for perm in perms:
        if perm.inv == 0:
            continue

        r = RCGraphRing()
        the_schub = ASx(perm, len(perm.trimcode)).change_basis(WordBasis)
        elem = 0
        for word, coeff in the_schub.items():
            cop = FA(*word).coproduct()
            for (cop1, cop2), coeff2 in cop.items():
                set1 = r.monomial(*cop1)
                set2 = r.monomial(*cop2)
                for rc1, rc2 in itertools.product(set1.keys(), set2.keys()):
                    elem += coeff * coeff2 * ASx(rc1.perm, len(rc1)) @ ASx(rc2.perm, len(rc2))

        print(elem)
            

        # the_front_schub = Schub(perm, len(perm.trimcode)).change_basis(MonomialBasis(PA.genset))
        # elem = 0
        # the_front_elem = 0
        # for perm2 in perms:
        #     for rc in RCGraph.all_rc_graphs(perm2, len(perm.trimcode)):
        #         # coeff = PA(*rc.length_vector).change_basis(SchubertPolyBasis(PA.genset)).get((perm, len(perm.trimcode)),0)
                
        #         # coeff = the_schub.get(rc.length_vector, 0)
        #         # if coeff != 0 and not rc.is_principal:
        #         #     print("Yay")
        #         elem += FA(*rc.length_vector) * the_schub.get(rc.length_vector, 0)
        #     #the_front_elem += PA(*rc.length_vector)# * FA(*rc.length_vector).change_basis(SchubertBasis).get((rc.perm, len(rc)), 0) * PA(*rc.length_vector).change_basis(SchubertPolyBasis(PA.genset)).get((rc.perm, len(rc)), 0)
        #     #w_coprod = FA(*rc.length_vector).coproduct()

        # assert elem == the_schub, f"Failed for {perm}, got elem={ elem} and the_schub={the_schub}"
        # #assert the_front_elem.almosteq(the_front_schub), f"Failed for {perm}, got front_elem={the_front_elem} and the_front_schub={the_front_schub}"
        # sys.exit(0)
        
        # print("schub")
        # elem_forest = 0
        # rc_prin = RCGraph.principal_rc(perm)
        # the_forest = ForestDual(*rc_prin.forest_weight).change_basis(SchubertBasis)
        
        # for theschub2, coeff0 in the_forest.items():
         
        #     for rc in RCGraph.all_rc_graphs(*theschub2):
        #         #if rc.forest_weight == rc_prin.forest_weight:
        #         coeff = ASx(*theschub2).change_basis(WordBasis).get(rc.length_vector, 0)
        #         #elem_forest += ForestDual(*rc.forest_weight) * PA(*rc.length_vector).change_basis(ForestPolyBasis(PA.genset)).get(rc.forest_weight,0)
        #         elem_forest += FA(*rc.length_vector) * coeff * coeff0
        # assert elem_forest == the_forest.change_basis(WordBasis), f"Failed for {perm}, got {elem_forest} and {the_forest}"
        # print("forest")