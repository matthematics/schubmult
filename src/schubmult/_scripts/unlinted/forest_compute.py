from functools import cache
from schubmult import *
from schubmult.rings.polynomial_algebra import *
from schubmult.symbolic import *
from schubmult.utils.perm_utils import weak_compositions
import itertools


if __name__ == "__main__":
    import sys

    n = int(sys.argv[1])

    if False:
        perms = Permutation.all_permutations(n)

        for perm in perms:
            if perm.inv == 0:
                continue
            polystink = {}
            for rc in RCGraph.all_rc_graphs(perm):
                inv = rc.forest_invariant
                polystink[inv] = polystink.get(inv, S.Zero) + rc.polyvalue(Sx.genset)
            for inv in polystink:
                # print(expand(polystink[inv] - Forest(*inv.forest.code).expand()) )
                # print(inv)
                # print(inv.forest.code)
                # print(Forest(*inv.forest.code).expand())
                assert expand(polystink[inv] - Forest(*inv.forest.code, 0).expand()) == 0
    from schubmult.rings.free_algebra import *
    ForestDual = FreeAlgebra(ForestBasis)
    r = RCGraphRing()
    invariant_stinkbat = {}
    product = {}
    rc_product = {}
    
    @cache
    def canonical_rc_rep(comp, basis):
        Farp = FreeAlgebra(basis)
        dualelem = Farp(*comp)
        signed_elem = r.from_free_algebra_element(dualelem)
        new_signed_elem = r.from_dict({**signed_elem})
        for rc, coeff in signed_elem.items():
            if coeff <= 0:
                new_signed_elem -= coeff * r(rc)
                next_rc = next(iter(rc0 for rc0, coeff2 in new_signed_elem.items() if rc0 != rc and rc0.perm == rc.perm and coeff == -coeff2))
                del new_signed_elem[next_rc]
            # if (rc.forest_invariant, len(rc)) not in invariant_stinkbat:
            #     invariant_stinkbat[(rc.forest_invariant, len(rc))] = rc
            # elif invariant_stinkbat[(rc.forest_invariant, len(rc))] != rc:
            #     continue
            # new_signed_elem[invariant_stinkbat[(rc.forest_invariant, len(rc))]] = new_signed_elem.get(invariant_stinkbat[(rc.forest_invariant, len(rc))], 0) + coeff
            # if invariant_stinkbat[(rc.forest_invariant, len(rc))] != rc:
            #     new_signed_elem[rc] = new_signed_elem.get(rc, 0) - coeff
        return new_signed_elem

    # for perm in perms:
    #     for length in range(len(perm.trimcode), n):
    #         for rc in RCGraph.all_rc_graphs(perm, length):
    #             if (rc.forest_invariant, len(rc)) not in invariant_stinkbat:
    #                 invariant_stinkbat[(rc.forest_invariant, len(rc))] = rc
    #             elif invariant_stinkbat[(rc.forest_invariant, len(rc))] != rc:
    #                 continue
    #             # coproduct = {}
    #             for cut in range(1, len(rc)):
    #                 rc1, rc2 = rc.vertical_cut(cut)
    #                 if (rc1.forest_invariant, len(rc1)) not in invariant_stinkbat:
    #                     invariant_stinkbat[(rc1.forest_invariant, len(rc1))] = rc1
    #                 elif invariant_stinkbat[(rc1.forest_invariant, len(rc1))] != rc1:
    #                     continue
    #                 if (rc2.forest_invariant, len(rc2)) not in invariant_stinkbat:
    #                     invariant_stinkbat[(rc2.forest_invariant, len(rc2))] = rc2
    #                 elif invariant_stinkbat[(rc2.forest_invariant, len(rc2))] != rc2:
    #                     continue
    #             #product[(rc1.forest_weight, rc2.forest_weight)] = coproduct.get((rc1.forest_weight, rc2.forest_weight), S.Zero) + ForestDual(*rc.forest_weight)
    #             #prd = Forest(*rc1.forest_weight) * Forest(*rc2.forest_weight)
    #             #rc_product[(rc1.forest_weight, rc2.forest_weight)] = rc_product.get((rc1.forest_weight, rc2.forest_weight), r.zero) + r(rc1) * r(rc2)
    #             rc_product[(rc1.forest_weight, rc2.forest_weight)] = r(rc1) * r(rc2)
    #             #assert rc_product.get(rc, 0) == prd.get(rc.forest_weight, 0), f"Failed for {perm} with cut {cut}, got {rc_product.get(rc, 0)}, expected {prd.get(rc.forest_weight, 0)}"
    # for (comp1, comp2), reduced_rc_elem in rc_product.items():
    #     forest_elem = sum([coeff * ForestDual(*rc.forest_weight) for rc, coeff in reduced_rc_elem.items()])
    #     assert forest_elem.almosteq(ForestDual(*comp1) * ForestDual(*comp2)), f"Failed for {comp1} and {comp2}, got {forest_elem}, expected {ForestDual(*comp1) * ForestDual(*comp2)}"
    from sympy import pretty_print
    compositions = weak_compositions(n, n)
    for comp in compositions:
        rc_elem = canonical_rc_rep(comp, KeyBasis)
        pretty_print(rc_elem)
        # leading_perm = next(iter(rc for rc in rc_elem if rc.is_principal)).perm
        # try_result = sum([coeff * PA(*rc.length_vector) for rc, coeff in rc_elem.items() if rc.perm == leading_perm])
        # assert try_result.almosteq(Forest(*comp).change_basis(MonomialBasis)), f"Failed for {comp}"
        # for comp2 in compositions:
        #     to_polify = Schub(uncode(comp2), len(comp2)).change_basis(MonomialBasis)
        #     assert rc_elem.to_free_algebra_element(WordBasis).pairing(to_polify) >= 0, f"Failed for {comp} and {comp2}, got {rc_elem.to_free_algebra_element().pairing(to_polify)}, expected nonnegative"


    print("Pinto bean exit")
    sys.exit()
    for comp1, comp2 in itertools.product(compositions, repeat=2):
        forest_elem1 = ForestDual(*comp1)
        forest_elem2 = ForestDual(*comp2)
        prod = forest_elem1 * forest_elem2
        rc_elem1 = canonical_rc_rep(comp1)
        rc_elem2 = canonical_rc_rep(comp2)
        rc_prod = rc_elem1 * rc_elem2
        assert prod.almosteq(rc_prod.to_free_algebra_element(ForestBasis)), f"Failed for {comp1} and {comp2}, got {rc_prod.to_free_algebra_element()}, expected {prod}"
        #assert rc_elem1.get(invariant_stinkbat[comp1], 0) == rc_elem2, f"Failed for {comp1} and {comp2}, got {rc_elem2}, expected {rc_elem1.get(invariant_stinkbat[comp1], 0)}"
    # for perm1, perm2 in itertools.product(perms, repeat=2):
    #     for length1, length2 in itertools.product(range(len(perm1.trimcode), n), range(len(perm2.trimcode), n)):
    #         for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1,length1), RCGraph.all_rc_graphs(perm2,length2)):
    #             inv1 = rc1.forest_invariant
    #             inv2 = rc2.forest_invariant
    #             poly1 = rc1.polyvalue(Sx.genset)
    #             poly2 = rc2.polyvalue(Sx.genset)
    #             prod = poly1 * poly2
    #             for k, v in prod.as_polynomial().items():
    #                 assert v >= 0, f"Negative coefficient {v} for monomial {k} in product of {perm1} and {perm2} at lengths {length1} and {length2}"