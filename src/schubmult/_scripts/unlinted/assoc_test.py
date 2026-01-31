from sympy import init_printing, pretty_print
from schubmult import *
from schubmult.rings.alt_rc_graph_ring import AltRCGraphRing

def is_isomorphic_crystal(tensor, rc_graph):
    hw1 = tensor.to_highest_weight()[0]
    hw2 = rc_graph.to_highest_weight()[0]
    if hw1.crystal_weight != hw2.crystal_weight:
        return False
    lw1 = hw1.to_lowest_weight()[0]
    lw2 = hw2.to_lowest_weight()[0]
    if lw1.crystal_weight != lw2.crystal_weight:
        return False
    return True


if __name__ == "__main__":
    # test module functionality

    import itertools
    import sys

    
    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)
    # sKEW DIV DIFF WEIGHT
    # is dual pieri Cauchy?
    r = RCGraphRing()
    for perm1, perm2, perm3 in itertools.product(perms, repeat=3):
        for length in range(1, n):
            length1 = length
            length2 = length
            length3 = length
            for rc1, rc2, rc3 in itertools.product(
                RCGraph.all_rc_graphs(perm1, length1),
                RCGraph.all_rc_graphs(perm2, length2),
                RCGraph.all_rc_graphs(perm3, length3),
            ):
                elem1 = ((r(rc1.transpose().normalize()) * r(rc2.transpose().normalize())) * r(rc3.transpose().normalize())).transpose(length)
                elem2 = (r(rc1.transpose().normalize()) * (r(rc2.transpose().normalize()) * r(rc3.transpose().normalize()))).transpose(length)
                elem1 = r.from_dict({k: v for k, v in elem1.items() if k.perm.inv == rc1.perm.inv + rc2.perm.inv + rc3.perm.inv})
                elem2 = r.from_dict({k: v for k, v in elem2.items() if k.perm.inv == rc1.perm.inv + rc2.perm.inv + rc3.perm.inv})
                assert all(v == 0 for v in (elem1 - elem2).values()), f"Failed associativity for {rc1}, {rc2}, {rc3}:\n{elem1}\n{elem2}"
                #print(f"Passed associativity for {rc1}, {rc2}, {rc3}")
                print(elem1)
                print(elem2)