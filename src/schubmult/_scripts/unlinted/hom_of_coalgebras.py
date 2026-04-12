from schubmult import *

r = RCGraphRing()
g = BoundedRCFactorAlgebra()

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
    
    for perm in perms:
        for size in range(len(perm.trimcode), n):
            for rc, cem_dict in RCGraph.full_CEM(perm, size, partition=tuple((~(perm.mul_dominant())).trimcode)).items():
                if rc.perm != perm:
                    continue
                tensor = g.from_tensor_dict(cem_dict, size).coproduct()
                rc_elem = _make_tensor_rc_element(tensor)
                rc_elem2 = r(rc).vertical_coproduct()
                assert rc_elem.almosteq(rc_elem2), f"Failure for {perm} with CEM {rc} and tensor {tensor}\nGot {rc_elem}, expected {rc_elem2}"
            print("Success for", perm, "with size", size)

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