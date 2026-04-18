from schubmult import *
from sympy import pretty_print

r = RCGraphRing()
g = BoundedRCFactorAlgebra()

def hom_rc_to_factor(rc):
    from schubmult.rings.polynomial_algebra import Schub, ElemSymPolyBasis

    elem_elem = ASx(rc.perm, len(rc)).change_basis(ElementaryBasis)
    #Schub(rc.perm, len(rc)).change_basis(ElemSymPolyBasis)
    # ASx(rc.perm, len(rc)).change_basis(ElementaryBasis)
    ret_elem = g.zero
    for (elem_tup, _), coeff in elem_elem.items():
        partition = tuple(reversed([index + 1 for index, a in enumerate(elem_tup) if a != 0]))
        ret_elem += coeff * g.from_tensor_dict(RCGraph.full_CEM(rc.perm, len(rc), partition=partition).get(rc, {}), size=len(rc.perm))
    return ret_elem

def _make_tensor_rc_element(tensor):
    result = (r@r).zero
    for (key1, key2), coeff in tensor.items():
        result += coeff * g(key1).to_rc_graph_ring_element() @ g(key2).to_rc_graph_ring_element()
    return result

if __name__ == "__main__":
    import sys
    import itertools
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    
    for perm1, perm2 in itertools.product(perms, repeat=2):
        for size1, size2 in itertools.product(range(len(perm1.trimcode), n), range(len(perm2.trimcode), n)):
            #for (rc1, cem_dict), (rc2, cem_dict2) in itertools.product(RCGraph.full_CEM(perm1, size1).items(), RCGraph.full_CEM(perm2, size2).items()):
            for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1, size1), RCGraph.all_rc_graphs(perm2, size2)):
                
                elem1 = hom_rc_to_factor(rc1)
                elem2 = hom_rc_to_factor(rc2)
                elem_prod = g.zero
                rc_prod =r(rc1) * r(rc2)
                for the_rc, coeff in rc_prod.items():
                    elem_prod += coeff * hom_rc_to_factor(the_rc)
                pretty_print(elem1)
                pretty_print(elem2)
                pretty_print(elem_prod)
                # rc_elem = _make_tensor_rc_element(tensor)
                # rc_elem2 = r(rc1).vertical_coproduct()
                # assert rc_elem.almosteq(rc_elem2), f"Failure for {perm1} with CEM {rc} and tensor {tensor}\nGot {rc_elem}, expected {rc_elem2}"
            print("Success for", perm1, "with size", size1)

    # for perm1, perm2 in itertools.product(perms, repeat=2):
    #     if perm1.inv == 0 or perm2.inv == 0:
    #         continue
    #     for size1, size2 in itertools.product(range(len(perm1.trimcode), n), range(len(perm2.trimcode), n)):
    #         for rc1, cem_dict1 in RCGraph.full_CEM(perm1, size1, partition=tuple((~(perm1.mul_dominant())).trimcode)).items():
    #             if rc1.perm != perm1:
    #                 continue
    #             for rc2, cem_dict2 in RCGraph.full_CEM(perm2, size2, partition=tuple((~(perm2.mul_dominant())).trimcode)).items():
    #                 if rc2.perm != perm2:
    #                     continue
    #                 tensor1 = next(iter(g.from_tensor_dict(cem_dict1, size1)))
    #                 tensor2 = next(iter(g.from_tensor_dict(cem_dict2, size2)))
    #                 tensor_prod = g.dual_product_on_basis(tensor1, tensor2)
    #                 rc_elem = tensor_prod.to_rc_graph_ring_element()
    #                 rc_elem2 = r(rc1)*r(rc2)
    #                 assert rc_elem.almosteq(rc_elem2), f"Failure for {perm1}, {perm2} with CEMs {rc1}, {rc2} and tensors {tensor1}, {tensor2}\nGot {rc_elem}, expected {rc_elem2}"
    #         print("Success for", perm1, perm2, "with sizes", size1, size2)