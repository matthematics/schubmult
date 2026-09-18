import itertools

from schubmult import CrystalGraph, CrystalGraphTensor
from schubmult import RCGraph
from schubmult import RCGraphRing

from schubmult import Permutation

# contactenating with the product gives new crystal graphs
# the tensor product is contained in it

if __name__ == "__main__":
    import sys

    from sympy import pretty_print
    rc_ring = RCGraphRing()
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(int(sys.argv[1]))
    for perm1, perm2 in itertools.product(perms, perms):
        for rc in RCGraph.all_rc_graphs(perm1):
            for rc2 in RCGraph.all_rc_graphs(perm2):
                rc_tensor = CrystalGraphTensor(rc.resize(n-1), rc2.resize(n-1))
                rc_elem_tensors = rc_ring.potential_products(rc2.resize(n-1), rc.resize(n-1), n-1)
                rc_tensor_set = set(rc_elem_tensors)

                pretty_print(rc_tensor)
                for i in range(1, n-1):
                    got_one_raise = False
                    got_one_lower = False
                    if rc_tensor.raising_operator(i):
                        for rct in rc_tensor_set:
                            if rct.length_vector != tuple([a+b for a, b in zip(rc.resize(n-1).length_vector, rc2.resize(n-1).length_vector)]):
                                continue
                            if rct.raising_operator(i) is not None:
                                got_one_raise = True
                        if not got_one_raise:
                            print(f"Missing raise at index {i} for tensor {rc_tensor} from {rc} and {rc2}")
                    if rc_tensor.lowering_operator(i):
                        for rct in rc_tensor_set:
                            if rct.length_vector != tuple([a+b for a, b in zip(rc.resize(n-1).length_vector, rc2.resize(n-1).length_vector)]):
                                continue
                            if rct.lowering_operator(i) is not None:
                                got_one_lower = True
                        if not got_one_lower:
                            print(f"Missing lower at index {i} for tensor {rc_tensor} from {rc} and {rc2}")
