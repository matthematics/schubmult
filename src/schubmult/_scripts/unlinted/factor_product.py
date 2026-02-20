from schubmult import *
from sympy import pretty_print, S
from schubmult.utils.schub_lib import elem_sym_perms_op, elem_sym_perms
from schubmult.combinatorics.crystal_graph import CrystalGraphTensor

def the_prod(tup):
    if len(tup) == 0:
        return RCGraph()
    return the_prod(tup[:-1]).resize(len(tup[-1])).squash_product(tup[-1])

if __name__ == "__main__":
    import sys
    import itertools

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    r = RCGraphRing()
    factorizations = {RCGraph(): ()}
    for k in range(1, n):
        new_factorizations = {}
        for p in range(0, k + 1):
            for rc, factor_tuple in factorizations.items():
                rc = rc.resize(k)
                for elem_rc in RCGraph.all_rc_graphs(uncode([0] * (k - p) + [1] * p), k):
                    new_rc = rc.squash_product(elem_rc)
                    new_factorizations[new_rc] = (*factor_tuple, elem_rc)
        factorizations = new_factorizations

    # divdiff_factorizations = {RCGraph(): {}}

    # for k in range(1, n):
    #     new_factorizations = {}
    #     for p in range(0, k + 1):
    #         for rc, factor_tuple in factorizations.items():
    #             rc = rc.resize(k)
    #             spot_perm =  (~(rc.perm)) * Permutation.w0(k)
    #             upwards = elem_sym_perms(spot_perm, k, k)
    #             for up_perm, diff in upwards:
    #                 if diff == 0:
    #                     new_rc = rc
    #                     new_factorizations[new_rc] = (*factor_tuple, RCGraph([]).resize(k))
    #                     continue
    #                 vec = [0] * k
    #                 spots = [i for i in range(len(up_perm)) if up_perm[i] == spot_perm[i] and i < k]
    #                 for j in spots:
    #                     vec[j] = 1
    #                 elem_rc = next(iter(RCGraph.all_rc_graphs(uncode([0] * (diff) + [1] * (k - diff)), k, weight=tuple(vec))))
    #                 new_rc = rc.squash_product(elem_rc)
    #                 test_perm = (~(new_rc.perm)) * Permutation.w0(k + 1)
    #                 assert test_perm == up_perm, f"Failed on {rc}\n{k=}\n{elem_rc}\n{new_rc}\n{up_perm=}\n{test_perm=}"
    #                 new_factorizations[new_rc] = (*factor_tuple, elem_rc)
    #     factorizations = new_factorizations
    # success = 0
    # for rc1, rc2 in itertools.product(factorizations.keys(), repeat=2):
    #     factor1 = factorizations[rc1]
    #     factor2 = factorizations[rc2]
    #     assert the_prod(factor1) == rc1, f"Failed on {rc1} with factors {factor1}"
    #     assert the_prod(factor2) == rc2, f"Failed on {rc2} with factors {factor2}"
    #     product_graph = RCGraph()
    #     edge = 0
    #     index = 0
    #     working_tup
    #     while index < len(factor2):
    #         while len(the_prod(factor1[:edge]).perm.trimcode)<=index + 1:
    #             edge += 1
    #     edge -= 1

    #     assert product_graph.perm in Sx(rc1.perm) * Sx(rc2.perm), f"Failed on \n{rc1} and \n{rc2}\n{product_graph.perm} not in {Sx(rc1.perm) * Sx(rc2.perm)}\n{product_graph}"
    #     success += 1
    # print(f"Successes: {success}")
    for perm in perms:
        for rc in RCGraph.all_rc_graphs(perm, n - 1):
            factor_tuple = factorizations[rc]
            tensor = factor_tuple[0]
            for i in range(1, len(factor_tuple)):
                tensor = CrystalGraphTensor(tensor, factor_tuple[i].resize(n-1))
            tensor, _ = tensor.to_lowest_weight()
            rc_lw = rc.to_lowest_weight()[0]
            assert rc_lw.crystal_weight == tensor.crystal_weight, f"Failed on {rc}\n{tensor}\n{rc_lw}"
            to_tuple = []
            the_factor = tensor
            for i in range(len(factor_tuple) - 1):
                to_tuple.insert(0, the_factor.factors[1])
                the_factor = the_factor.factors[0]
            to_tuple.insert(0, the_factor)
            product = RCGraph()
            for i, rct in enumerate(to_tuple):
                rct = rct.resize(i + 1)
                product = product.resize(len(rct)).squash_product(rct)
            assert product == rc.to_lowest_weight()[0], f"Failed on {rc}\n{to_tuple}\n{product}"

            tensor, _ = tensor.to_highest_weight()

            to_tuple = []
            the_factor = tensor
            for i in range(len(factor_tuple) - 1):
                to_tuple.insert(0, the_factor.factors[1])
                the_factor = the_factor.factors[0]
            to_tuple.insert(0, the_factor)
            product = RCGraph()
            for i, rct in enumerate(to_tuple):
                rct = rct.resize(i + 1)
                product = product.resize(len(rct)).squash_product(rct)
            assert product == rc.to_highest_weight()[0], f"Failed on {rc}\n{to_tuple}\n{product}"
            # assert product.to_highest_weight()[0] == rc.to_highest_weight()[0], f"Failed on {rc}\n{to_tuple}\n{product}"
                
            