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

def convert_tensor(tensor):
        ret = (ASx@FA@ASx@FA).zero
        for (rc1, rc2), coeff in tensor.items():
            ret += coeff * ASx(rc1.perm, len(rc1)) @ FA(*rc1.length_vector) @ ASx(rc2.perm, len(rc2)) @ FA(*rc2.length_vector)
        return ret

def schub_weight_equal(tensor1, tensor2):
    
    
    return convert_tensor(tensor1).almosteq(convert_tensor(tensor2))


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

def full_coprod(perm, weight, length):
    cem = RCGraph.full_CEM(perm, length)
    result = 0
    for rc, cem_dict in cem.items():
        if rc.length_vector != weight:
            continue
        for rc_tup, coeff in cem_dict.items():
            if coeff == 0:
                continue
            if tuple([a + b for a, b in zip(rc_tup[0].length_vector, rc_tup[1].length_vector)]) == weight:
                yield (rc_tup, coeff)

if __name__ == "__main__":
    from schubmult.combinatorics.crystal_graph import CrystalGraphTensor
    from schubmult.abc import *
    from schubmult.rings.polynomial_algebra import *
    from sympy import pretty_print
    n = 5
    
    perms = Permutation.all_permutations(n)
    ring = GrassTensorAlgebra()
    

    grass_tensor_elems = {}
    for perm in perms:
        grass_tensor_elems[perm] = ring.zero
        cem = RCGraph.full_CEM(perm, n - 1)
        for rc, cem_dict in cem.items():
            grass_tensor_elems[perm] += ring.from_dict(cem_dict)
    for perm1, perm2 in itertools.product(perms, repeat=2):
        prd = Sx(perm1) * Sx(perm2)
        prd_elem = grass_tensor_elems[perm1] * grass_tensor_elems[perm2]
        prd_rc = prd_elem.to_rc_graph_ring_element()
        for rc, coeff in prd_rc.items():
            assert prd.get(rc.perm, 0) == coeff, f"Failed on {perm1} * {perm2}, got {rc.perm}: {coeff} which is not in {prd}\n{prd.get(rc.perm,0)=}\n{rc.perm=}"
        print(f"Checked {perm1} * {perm2}")