from schubmult import *

if __name__ == "__main__":
    from schubmult.combinatorics.crystal_graph import CrystalGraphTensor
    from schubmult.rings.polynomial_algebra import *
    import sys
    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)

    for perm in perms:
        for k in range(1, n):
            elem_rc = RCGraph.principal_rc(~uncode([k]), n - 1)
            the_dubya = Sx(perm) * Sx(elem_rc.perm)
            rc_set = set()
            for w in the_dubya.keys():
                rc_set.update(RCGraph.all_hw_rcs(w, n - 1))
            hw_map = {}
            prin_rc = RCGraph.principal_rc(perm, n - 1)
            omega_class = prin_rc.forest_invariant
            omega_decomp = {}
            for rc in RCGraph.all_rc_graphs(perm, n - 1):
                if rc.forest_invariant != omega_class:
                    continue
                tensor = CrystalGraphTensor(elem_rc, rc)
                # extremal_weight = tuple([a + b for a, b in zip(rc.extremal_weight, elem_rc.length_vector)])
                hw_tensor, raise_seq = tensor.to_highest_weight()
                if hw_tensor not in hw_map:
                    def dom_key(weight):
                        return tuple([sum(weight[:i]) for i in range(1, len(weight) + 1)])
                    extremal_weight = min([tens_bag.crystal_weight for tens_bag in hw_tensor.full_crystal if tuple(hw_tensor.crystal_weight) == tuple(sorted(tens_bag.crystal_weight, reverse=True))], key=lambda x: dom_key(x))
                    me_map = next(iter([rc2 for rc2 in rc_set if tuple(rc2.extremal_weight) == tuple(extremal_weight)]))
                    rc_set.remove(me_map)
                    hw_map[hw_tensor] = me_map
                mapped_rc = hw_map[hw_tensor].reverse_raise_seq(raise_seq)
                omega_decomp[mapped_rc.forest_invariant] = omega_decomp.get(mapped_rc.forest_invariant, set())
                omega_decomp[mapped_rc.forest_invariant].add(mapped_rc)
            for invariant, rcs in omega_decomp.items():
                rc_sum = sum([rc.polyvalue(Sx.genset) for rc in rcs])
                poly = PA.from_expr(rc_sum, length=n-1).change_basis(ForestPolyBasis)
                assert all(v >= 0 for v in poly.values()), f"Error: Invariant {invariant} has negative coefficient in decomposition: {poly}, {rcs=}"
                assert len([key for key, v in poly.items() if v > 0]) == 1, f"Error: Invariant {invariant} has more than one term in decomposition: {poly}"
                if any(v > 1 for v in poly.values()):
                    print("Sweesh it's notnin")
            # assert len(rc_set) == 0, f"Error: Not all RC graphs for {perm} are accounted for: missing {rc_set}"