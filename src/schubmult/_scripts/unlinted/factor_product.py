from schubmult import *
from sympy import pretty_print, S

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
    success = 0
    for rc1, rc2 in itertools.product(factorizations.keys(), repeat=2):
        factor1 = factorizations[rc1]
        factor2 = factorizations[rc2]
        assert the_prod(factor1) == rc1, f"Failed on {rc1} with factors {factor1}"
        assert the_prod(factor2) == rc2, f"Failed on {rc2} with factors {factor2}"
        product_graph = RCGraph()
        for i in range(len(factor1)):
            product_graph = product_graph.resize(i+1)
            product_graph = product_graph.squash_product(factor1[i]).squash_product(factor2[i])
        assert product_graph.perm in Sx(rc1.perm) * Sx(rc2.perm), f"Failed on \n{rc1} and \n{rc2}\n{product_graph.perm} not in {Sx(rc1.perm) * Sx(rc2.perm)}\n{product_graph}"
        success += 1
    print(f"Successes: {success}")