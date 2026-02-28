from schubmult import *
from schubmult.rings.free_algebra import *
from schubmult.utils.tuple_utils import *
from itertools import product
from sympy import pretty_print

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)
    seqs = [pad_tuple(perm.trimcode, n-1) for perm in perms]

    Forest = FreeAlgebra(ForestBasis)
    Key = FreeAlgebra(KeyBasis)
    r = RCGraphRing()

    rc_by_invariant = {}
    for perm in perms:
        rc_by_invariant[perm] = rc_by_invariant.get(perm, {})
        for rc in RCGraph.all_rc_graphs(perm):
            rc_by_invariant[perm][rc.forest_invariant] = rc_by_invariant[perm].get(rc.forest_invariant, set())
            rc_by_invariant[perm][rc.forest_invariant].add(rc)
    
    for perm1, perm2 in product(perms, repeat=2):
        for (invariant1, rc_set1), (invariant2, rc_set2) in product(rc_by_invariant[perm1].items(), rc_by_invariant[perm2].items()):
            # if rc1.perm in perms_seen and perms_seen[rc1.perm] != rc1:
            #     continue
            # if rc2.perm in perms_seen and perms_seen[rc2.perm] != rc2:
            #     continue
            # perms_seen[rc1.perm] = rc1
            # perms_seen[rc2.perm] = rc2
            rc1 = next(iter(rc_set1))
            rc2 = next(iter(rc_set2))
            forest1 = Forest(*rc1.forest_weight)
            forest2 = Forest(*rc2.forest_weight)
            r_elem1 = r(rc1)
            r_elem2 = r(rc2)
            # assert r_elem1.to_free_algebra_element(ForestBasis).almosteq(forest1)
            # assert r_elem2.to_free_algebra_element(ForestBasis).almosteq(forest2), f"Error: RC graph element for invariant {invariant1} or {invariant2} does not match its corresponding forest basis element for permutations {perm1} and {perm2}: expected {forest1} and {forest2}, got {r_elem1.to_free_algebra_element(ForestBasis).change_basis(ForestBasis)} and {r_elem2.to_free_algebra_element(ForestBasis).change_basis(ForestBasis)}"
            f_prod = forest1 * forest2
            r_prod = r_elem1 * r_elem2
            new_forest = 0
            invariants_seen = set()
            for rc3, coeff in r_prod.items():
                # if rc3.forest_weight != rc3.length_vector:
                #     continue
                if rc3.forest_invariant in invariants_seen:
                    continue
                invariants_seen.add(rc3.forest_invariant)
                # if rc3.forest_weight in invariants_seen_weight and invariants_seen_weight[rc3.forest_weight] != rc3.forest_invariant:
                #     continue
                # invariants_seen_weight[rc3.forest_weight] = rc3.forest_invariant
                assert coeff == 1
                #new_forest += coeff * ASx(rc3.perm, len(rc3)).change_basis(ForestBasis)
                new_forest += coeff * Forest(*rc3.forest_weight)
            assert f_prod == new_forest, f"Error: Forest product does not match RC product for {rc1} and {rc2}: {f_prod} != {new_forest}\nRC product: {r_prod}"
    #     # rc1 = RCGraph.principal_rc(perm1)
    #     # rc2 = RCGraph.principal_rc(perm2)
    #         # if rc1.length_vector != pad_tuple(rc1.forest_invariant.forest.code, len(rc1)) or rc2.length_vector != pad_tuple(rc2.forest_invariant.forest.code, len(rc2)):
    #         #     #print(f"Warning: RC graph {rc1} or {rc2} has length vector that does not match its forest invariant. Skipping.")
    #         #     continue
    #         forest1 = Forest(*pad_tuple(rc1.forest_invariant.forest.code, len(rc1)))
    #         forest2 = Forest(*pad_tuple(rc2.forest_invariant.forest.code, len(rc2)))
    #         r_elem1 = r(rc1)
    #         r_elem2 = r(rc2)
    #         f_prod = forest1 * forest2
    #         r_prod = r_elem1 * r_elem2
    #         new_forest = 0
    #         invars_seen = set()
    #         #perms_seen = set()
    #         for rc3, coeff in r_prod.items():
    #             # if pad_tuple(rc3.forest_invariant.forest.code, len(rc3)) != rc3.length_vector:
    #             #     #print(f"Warning: RC graph {rc3} has forest invariant that does not match its principal RC graph's forest invariant. Skipping.")
    #             #     continue
    #             #invars_seen.add(rc3.forest_invariant)
    #             if rc3.forest_invariant.forest.code in invars_seen:
    #                 continue
    #             invars_seen.add(rc3.forest_invariant.forest.code)
    #             assert coeff == 1
    #             new_forest += coeff * Forest(pad_tuple(rc3.forest_invariant.forest.code, len(rc3)))
    #         assert f_prod == new_forest, f"Error: Forest product does not match RC product for {rc1} and {rc2}: {f_prod} != {new_forest}\nRC product: {r_prod}"
    # print("Forest")
    # # for seq1, seq2 in product(seqs, repeat=2):
    # #     print("Sequence:", seq1, seq2)
    # #     #forest_elem = FA(*seq).change_basis(ForestBasis)
    # #     m1 = r.monomial(*seq1)
    # #     m2 = r.monomial(*seq2)

    # #     for rc1, rc2 in product(m1.keys(), m2.keys()):
    # #         # if rc1.length_vector != pad_tuple(rc1.forest_invariant.forest.code, len(rc1)) or rc2.length_vector != pad_tuple(rc2.forest_invariant.forest.code, len(rc2)):
    # #         #     #print(f"Warning: RC graph {rc1} or {rc2} has length vector that does not match its forest invariant. Skipping.")
    # #         #     continue
    # #         r_elem1 = r(rc1)
    # #         r_elem2 = r(rc2)
    # #         f_elem1 = Forest(*pad_tuple(rc1.forest_invariant.forest.code, len(rc1)))
    # #         f_elem2 = Forest(*pad_tuple(rc2.forest_invariant.forest.code, len(rc2)))
    # #         r_prod = r_elem1 * r_elem2
    # #         f_prod = f_elem1 * f_elem2
    # #         for rc3, coeff in r_prod.items():
    # #             # if pad_tuple(rc3.forest_invariant.forest.code, len(rc3)) != rc3.length_vector:
    # #             #     #print(f"Warning: RC graph {rc3} has forest invariant that does not match its principal RC graph's forest invariant. Skipping.")
    # #             #     continue
    # #             f_coeff = f_prod[pad_tuple(rc3.forest_invariant.forest.code, len(rc3))]
    # #             assert coeff == f_coeff, f"Error: Forest product does not match RC product for {rc1} and {rc2}: {f_coeff} != {coeff}"
    
    # forest_basis = {}
    # from schubmult.rings.polynomial_algebra import Forest, MonomialBasis
    # from schubmult.symbolic import expand_seq
    # ForestDual = FreeAlgebra(ForestBasis)
    # for seq in seqs:
    #     print("Sequence:", seq)
    #     forest_elem = ForestDual(*seq)
    #     pretty_print(r.from_free_algebra_element(forest_elem))
        # r_elem = r.from_free_algebra_element(forest_elem)
        # r_elem2 = r.from_dict({**r_elem})
        # for rc, coeff in r_elem.items():
        #     if coeff < 0:
        #         sub_elem = min([k for k, v in r_elem2.items() if v > 0 and k.perm == rc.perm], key=lambda k: pad_tuple(k.forest_invariant.forest.code, len(k)))
        #         del r_elem2[sub_elem]
        #         del r_elem2[rc]
        # pretty_print(r_elem2)
        # forest_basis[seq] = r_elem2

        
        # forest_poly = Forest(*seq).change_basis(MonomialBasis)
        # for monom, coeff in forest_poly.items():
        #     schub_monom = Sx(expand_seq(monom, forest_poly.ring.genset))
        #     r_elem_monom = r.zero
        #     for perm2, coeff2 in schub_monom.items():
        #         r_elem_monom += coeff2 * r.schub(perm2, len(monom))
        #     result = 0
        #     for rc1, coeff1 in r_elem_monom.items():
        #         for rc2, coeff2 in r_elem2.items():
        #             if rc1.perm == rc2.perm:
        #                 result += coeff1 * coeff2
        #     assert result == coeff, f"Error: {result=} Forest basis element {r_elem_monom} does not have inner product {coeff} with its corresponding RC graph element: {forest_basis[seq]}"
        #     the_monomial
        #     if seq == seq2:
        #         assert forest_basis[seq].to_free_algebra_element().pairing(forest_poly) == 1, f"Error: Forest basis element {forest_elem} does not have inner product 1 with its corresponding RC graph element: {forest_basis[seq]}"
        #     else:
        #         assert forest_basis[seq].to_free_algebra_element().pairing(forest_poly) == 0, f"Error: Forest basis element {forest_elem} does not have inner product 0 with non-corresponding RC graph element {forest_poly}: {forest_basis[seq]}"
        #print("PAIRIARI")

    # for seq in seqs:
    #     print("Sequence:", seq)
    #     forest_elem = 
        
        # print(set((pad_tuple(rc.forest_invariant.forest.code, len(rc))) for rc in r_elem2.keys()))
        # print([(coeff,rc.forest_invariant) for rc, coeff in r_elem.items()])
    # from schubmult.rings.polynomial_algebra import Schub, PA, ForestPolyBasis, PolynomialAlgebra
    # ForestPoly = PolynomialAlgebra(ForestPolyBasis)
    # PA.pa
    # for seq in seqs:
    #     schub_elem = FA(*seq).change_basis(SchubertBasis)
    #     for forest_seq in seqs:
    #         coeff = 0
    #         for rc in forest_basis[forest_seq].keys():
    #             if (rc.perm, len(rc)) in schub_elem:
    #                 coeff += schub_elem[(rc.perm, len(rc))]
                
    #         forest_monom = schub_elem.poly_inner_product(forest_basis[forest_seq])

    # # print("Key")
    # # for seq in seqs:
    # #     print("Sequence:", seq)
    # #     r_elem = r.from_free_algebra_element(Key(*seq))
    # #     interm = {k: sum([v2 for k2, v2 in r_elem.items() if k2.perm == k.perm]) for k, v in r_elem.items()}
    # #     r_elem = r.from_dict({k: v for k, v in interm.items() if v != 0})
    # #     pretty_print(r_elem)
    # #     print([(coeff,rc.extremal_weight) for rc, coeff in r_elem.items()])
        