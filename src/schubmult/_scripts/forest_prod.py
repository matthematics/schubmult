from schubmult import *
from schubmult.rings.polynomial_algebra import *
from schubmult.rings.free_algebra import *
from schubmult.rings.combinatorial.forest_rc_ring import ForestRCGraphRing
from schubmult.combinatorics.indexed_forests import *
from sympy import pretty_print

# def snap_key(rc):
#     return next(iter([rcc for rcc in RCGraph.all_rc_graphs(rc.perm, len(rc), weight=rc.extremal_weight)]))

# def forest_equal(rc1, rc2):
#     return rc1.forest_weight == rc2.forest_weight and rc1.omega_invariant[1] == rc2.omega_invariant[1]

CompSchub = FreeAlgebra(CompositionSchubertBasis)
def omega_invariant(word):
    from schubmult.combinatorics.indexed_forests import letterpair, omega_insertion

    def word_to_pair_labeled(word):
        counts = {}
        out = []
        for a in word:
            aa = int(a)
            counts[aa] = counts.get(aa, 0) + 1
            out.append(letterpair(aa, counts[aa]))
        return tuple(out)

    return omega_insertion(word_to_pair_labeled(word))

def is_in_forest_ideal(rc_graph_ring_element):
    buildup_dict = {}
    for rc, coeff in rc_graph_ring_element.items():
        if coeff == 0:
            continue
        buildup_dict[rc.forest_weight] = buildup_dict.get(rc.forest_weight, 0) + coeff
    return len({k: v for k, v in buildup_dict.items() if v != 0}) == 0

r = RCGraphRing()
def forestify(rc_graph_ring_element):
    from schubmult.combinatorics.indexed_forests import IndexedForest

    out = 0
    for rc, coeff in rc_graph_ring_element.items():
        if coeff == 0:
            continue
        if rc.forest_weight == rc.perm.pad_code(len(rc)):
            #continue
        #if rc.perm.pad_code(len(rc)) == rc.forest_weight:
            out += coeff * ForestDual(*rc.forest_weight)#.change_basis(CompositionSchubertBasis).get(tuple(rc.perm.pad_code(len(rc))), 0) * CompSchub(*rc.perm.pad_code(len(rc))) #@ CompSchub(*rc.perm.pad_code(len(rc))) @ FA(*rc.length_vector)
        #if tuple(rc.perm.pad_code(len(rc))) in ForestDual(*rc.forest_weight).change_basis(CompositionSchubertBasis):
        # else:
        #     out += coeff * ForestDual(*rc.forest_weight)
        # else:
        #     out += coeff * ForestDual(*rc.perm.pad_code(len(rc))) #@ CompSchub(*rc.perm.pad_code(len(rc))) @ FA(*rc.length_vector)
        #CompSchub(rc.perm.pad_code(len(rc)))
        #ASx.from_dict({(k.perm, len(rc)): 1 for k, v in r.from_free_algebra_element(ForestDual(*rc.forest_weight)).items() if k.forest_weight == rc.forest_weight and k.perm == rc.perm})
    return out

if __name__ == "__main__":
    import sys
    import itertools

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    r = RCGraphRing()

    for perm1, perm2 in itertools.product(perms, repeat=2):
        for length1, length2 in itertools.product(range(max(1,len(perm1.trimcode)), n), range(max(1,len(perm2.trimcode)), n)):
            comp1 = perm1.pad_code(length1)
            comp2 = perm2.pad_code(length2)
            rep1 = ForestDual(*comp1).change_basis(SchubertBasis)
            rep2 = ForestDual(*comp2).change_basis(SchubertBasis)
            test_item = 0
            for p1, p2 in itertools.product(rep1.keys(), rep2.keys()):
                
                expected = 0
                for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(p1[0], p1[1]), RCGraph.all_rc_graphs(p2[0], p2[1])):
                    forest_rc1 = forestify(r(rc1))
                    forest_rc2 = forestify(r(rc2))
                    if forest_rc1 == 0 or forest_rc2 == 0:
                        continue
                    test_item += forestify(r(rc1) * r(rc2))
                    break
                #break
            expected = ForestDual(*comp1) * ForestDual(*comp2)#+= forestify(r(rc1)) * forestify(r(rc2))
                    #assert is_in_forest_ideal(test_item), f"Failed for {perm1} and {perm2} with lengths{length1} and {length2}, got {test_item}"
            assert test_item.almosteq(expected), f"Failed for {perm1} and {perm2} with lengths{length1} and {length2}, got {test_item} but expected {expected}\ndiff = {test_item - expected}"
            print(f"Success for {perm1} and {perm2} with lengths{length1} and {length2}")
        # by_weight = {}
        # for rc in RCGraph.all_rc_graphs(perm, n - 1):
        #     by_weight[rc.forest_weight] = by_weight.get(rc.forest_weight, set()) | {rc}
        # for weight, rc_set in by_weight.items():
        #     for rc1, rc2 in itertools.product(rc_set, repeat=2):
        #         test_item = r(rc1) - r(rc2)
        #         assert is_in_forest_ideal(test_item), f"Failed for {perm} with weight {weight}, got {test_item}"
        #         for perm2 in perms:
        #             for rc in RCGraph.all_rc_graphs(perm2):
        #                 test_item_left = r(rc) * test_item
        #                 assert is_in_forest_ideal(test_item_left), f"Not left ideal"
        #                 test_item_right = test_item * r(rc)
        #                 assert is_in_forest_ideal(test_item_right), f"Not right ideal"
        #     print(f"Success for {perm} with weight {weight}")
        