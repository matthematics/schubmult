from schubmult import *
from schubmult.rings.polynomial_algebra import *
from schubmult.rings.free_algebra import FreeAlgebra, ForestBasis, KeyBasis
from schubmult.utils.tuple_utils import pad_tuple
from schubmult.utils.perm_utils import add_perm_dict

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Verify RC-graph transitions for forest and key invariants.")
    parser.add_argument("n", type=int, help="Permutation size.")
    parser.add_argument("--skip-forest", action="store_true", help="Skip forest-basis verification.")
    parser.add_argument("--skip-key", action="store_true", help="Skip key-basis verification.")
    args = parser.parse_args()

    n = args.n
    perms = Permutation.all_permutations(n)
    if not args.skip_forest:
        FS = FreeAlgebra(ForestBasis)
        for perm in perms:
            forest_decomp = {}
            for rc in RCGraph.all_rc_graphs(perm):
                forest_decomp[rc.forest_invariant] = forest_decomp.get(rc.forest_invariant, set())
                forest_decomp[rc.forest_invariant].add(rc)
            for invariant, rcs in forest_decomp.items():
                forest_elem = PA.from_expr(sum([rc.polyvalue(PA.genset) for rc in rcs])).change_basis(Forest._basis)
                the_inv_elem = Forest(*invariant.forest.code)
                assert forest_elem.almosteq(the_inv_elem), f"Error: Invariant {invariant} does not match sum of RC graph polynomials for permutation {perm}: expected {the_inv_elem}, got {forest_elem}"

                for i in range(1, len(perm.trimcode) - 1):
                    coproduct = {}
                    invariant_to_rc = {}
                    bottom_invariant_to_rc = {}
                    for rc in rcs:
                        top_rc, bottom_rc = rc.vertical_cut(i)
                        top_invariant = top_rc.forest_invariant
                        bottom_invariant = bottom_rc.forest_invariant
                        if top_invariant not in invariant_to_rc:
                            invariant_to_rc[top_invariant] = top_rc
                        elif invariant_to_rc[top_invariant] != top_rc:
                            continue
                        if bottom_invariant not in bottom_invariant_to_rc:
                            bottom_invariant_to_rc[bottom_invariant] = bottom_rc
                        elif bottom_invariant_to_rc[bottom_invariant] != bottom_rc:
                            continue

                        coproduct[(pad_tuple(top_invariant.forest.code, i), pad_tuple(bottom_invariant.forest.code, len(perm.trimcode) - i))] = coproduct.get((pad_tuple(top_invariant.forest.code, i), pad_tuple(bottom_invariant.forest.code, len(perm.trimcode) - i)), 0) + 1

                    for (top_code, bottom_code), count in coproduct.items():
                        top_elem = FS(*top_code)
                        bottom_elem = FS(*bottom_code)
                        prod_elem = top_elem * bottom_elem
                        assert prod_elem.get(pad_tuple(invariant.forest.code, len(perm.trimcode)), 0) == count, f"Error: Invariant {invariant} has incorrect count in product of {top_code} and {bottom_code}: expected {count}, got {prod_elem.get(pad_tuple(invariant.forest.code, len(perm.trimcode)), 0)}, rcs={rcs}"
        print("Forest verified")

    if not args.skip_key:
        KS = FreeAlgebra(KeyBasis)
        Key = PolynomialAlgebra(KeyPolyBasis(Sx.genset))
        for perm in perms:
            key_decomp = {}
            for rc in RCGraph.all_rc_graphs(perm):
                invariant = rc.to_highest_weight()[0]
                key_decomp[invariant] = key_decomp.get(invariant, set())
                key_decomp[invariant].add(rc)
            for invariant, rcs in key_decomp.items():
                key_elem = PA.from_expr(sum([rc.polyvalue(PA.genset) for rc in rcs]), length=len(perm.trimcode)).change_basis(Key._basis)
                the_inv_elem = Key(*invariant.extremal_weight)
                assert key_elem.almosteq(the_inv_elem), f"Error: Invariant {invariant} does not match sum of RC graph polynomials for permutation {perm}: expected {the_inv_elem}, got {key_elem}"

                for i in range(1, len(perm.trimcode) - 1):
                    coproduct = {}
                    invariant_to_rc = {}
                    bottom_invariant_to_rc = {}
                    for rc in rcs:
                        if not rc.is_highest_weight:
                            continue
                        top_rc, bottom_rc = rc.vertical_cut(i)
                        top_invariant = top_rc.to_highest_weight()[0]
                        bottom_invariant = bottom_rc.to_highest_weight()[0]
                        if top_invariant not in invariant_to_rc:
                            invariant_to_rc[top_invariant] = top_rc
                        elif invariant_to_rc[top_invariant] != top_rc:
                            continue
                        if bottom_invariant not in bottom_invariant_to_rc:
                            bottom_invariant_to_rc[bottom_invariant] = bottom_rc
                        elif bottom_invariant_to_rc[bottom_invariant] != bottom_rc:
                            continue

                        coproduct[(top_invariant.extremal_weight, bottom_invariant.extremal_weight)] = coproduct.get((top_invariant.extremal_weight, bottom_invariant.extremal_weight), 0) + 1

                    for (top_code, bottom_code), count in coproduct.items():
                        top_elem = KS(*top_code)
                        bottom_elem = KS(*bottom_code)
                        prod_elem = top_elem * bottom_elem
                        assert prod_elem.get(invariant.extremal_weight, 0) == count, f"Error: Invariant {invariant} has incorrect count in product of {top_code} and {bottom_code}: expected {count}, got {prod_elem.get(invariant.extremal_weight, 0)}, rcs={rcs}"
        print("Key verified")

    if args.skip_forest and args.skip_key:
        print("Nothing to do: both --skip-forest and --skip-key were set.")
                
        