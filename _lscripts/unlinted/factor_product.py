from schubmult import *
from sympy import pretty_print, S
from schubmult.combinatorics.crystal_graph import CrystalGraphTensor
# from schubmult.symbolic import expand_seq
# from schubmult.rings.polynomial_algebra import *

# def the_prod(tup):
#     if len(tup) == 0:
#         return RCGraph()
#     return the_prod(tup[:-1]).resize(len(tup[-1])).squash_product(tup[-1])

if __name__ == "__main__":
    import sys
    import itertools

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    r = DualRCGraphRing()
    #factorizations = {RCGraph(): set([((),())])}
    for perm in perms:
        for rc in RCGraph.all_hw_rcs(perm, len(perm.trimcode)):
            # get factorization
            tensor = CrystalGraphTensor(*rc.squash_decomp())
            stack = [(tensor, rc)]
            while stack:
                tensor, rc = stack.pop()
                for i in range(1, len(rc)):
                    rc_check = rc.lowering_operator(i)
                    tensor_check = tensor.lowering_operator(i)
                    if rc_check is None:
                        assert tensor_check is None, f"Failure for {rc}, {tensor}, {i}"
                        continue
                    assert tensor_check is not None, f"Failure for {rc}, {tensor}, {i}"
                    stack.append((tensor_check, rc_check))
            print("rc succefefeis")
            