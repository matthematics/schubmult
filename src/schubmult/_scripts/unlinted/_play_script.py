from multiprocessing import Pool, cpu_count
from schubmult.rings.polynomial_algebra import *
from schubmult import *
from schubmult.abc import x
from schubmult.utils.perm_utils import mu_A
from schubmult.rings.free_algebra import *
from schubmult.utils.tuple_utils import pad_tuple
import itertools

def check_perm(perm_arr):
    perm = Permutation(perm_arr)
    n = len(perm_arr)
    r = RCGraphRing()
    for i in range(len(perm.trimcode), n):
        rc_elem = r.from_free_algebra_element(ASx(perm, i))
        if not all(int(v) in {1, 0, -1} for v in rc_elem.values()):
            return (False, perm_arr, dict(rc_elem))
    return (True, perm_arr, None)

EP = PolynomialAlgebra(ElemSymPolyBasis(Sx.genset))

# def in_good_algebra(elem_comp, monomial_vec):
#     #monom = EP.from_dict({elem_comp: 1}).change_basis(EP._basis.monomial_basis)
#     #return monom.get(monomial_vec, 0) != 0
#     assert isinstance(elem_comp[1], Permutation)
#     if FA(*)
#     return True

SS = FreeAlgebra(SchubertSchurBasis)
EE = FreeAlgebra(ElementaryBasis)


r = RCGraphRing()
def rc_to_ss(rc1):
    ret = EE.zero
    if isinstance(rc1, RCGraph):
        schuby = r(rc1).to_free_algebra_element().change_basis(ElementaryBasis)
        wordy = FA(*rc1.length_vector).change_basis(ElementaryBasis)
        
        for elem_comp, coeff in schuby.items():
            if coeff == 0:
                continue
            if wordy.get(elem_comp, 0) != 0:
                ret += coeff * EE(*elem_comp)
    else:
        for rc, coeff in rc1.items():
            ret += coeff * rc_to_ss(rc)
    return ret

if __name__ == "__main__":
    
    n = 4
    
    perms = Permutation.all_permutations(n)
    
    #print(Schub((uncode([0,2,1]),3)).change_basis(ElemSymPolyBasis))
    for perm1, perm2 in itertools.product(perms, repeat=2):
        result = 0
        
        prd = Sx(perm1) * Sx(perm2)
        #perms_seen = {}
        seen = set()
        for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1, n -1), RCGraph.all_rc_graphs(perm2, n - 1)):
            wt = tuple([a + b for a, b in zip(rc1.length_vector, rc2.length_vector)])
            if (rc1.length_vector, rc2.length_vector) in seen or (rc2.length_vector, rc1.length_vector) in seen:
                continue
            seen.add((rc1.length_vector, rc2.length_vector))
               
            result += Schub.from_dict(r.monomial(*wt).to_free_algebra_element())
                
                    #perms_seen[wt] = perms_seen.get(wt, 0)  + Sx(w)
            
        assert result == Schub(perm1, n - 1) * Schub(perm2, n - 1),  f"Failed on {perm1} * {perm2}, got {result}"
        # if perm.inv == 0:
        #     continue
        # for rc in RCGraph.all_rc_graphs(perm):
        #     free_schub_elem = r(rc).to_free_algebra_element().change_basis(SchubertSchurBasis)
        #     accum_free = 0
        #     for elem_comp, coeff in free_schub_elem.items():
        #         #print(elem_comp)s
        #         # elem_monoms = EP.from_dict({elem_comp: 1}).change_basis(EP._basis.monomial_basis)
        #         # if elem_monoms.get(rc.length_vector, 0) != 0:
        #         #     accum_free += coeff * EE(*elem_comp)
        #         if in_good_algebra(elem_comp, rc.length_vector):
        #             accum_free += SS(*elem_comp)

        #     tensor_elem = accum_free @ FA(*rc.length_vector)
        #     new_tensor_elem = 0
        #     for (elem_comp1, monom_comp1, elem_comp2, monom_comp2), coeff in tensor_elem.coproduct().items():
        #         if coeff != 0:
        #             # if in_good_algebra(elem_comp1, monom_comp1) and in_good_algebra(elem_comp2, monom_comp2):
        #             #     new_tensor_elem += coeff * (EE(*elem_comp1) @ FA(*monom_comp1)) * (EE(*elem_comp2) @ FA(*monom_comp2))

        #             if sum(elem_comp1[0]) == sum(monom_comp1):
        #                 assert in_good_algebra(elem_comp1, monom_comp1)

        #             if sum(elem_comp2[0]) == sum(monom_comp2):
        #                 assert in_good_algebra(elem_comp2, monom_comp2)
        #     print(new_tensor_elem)
