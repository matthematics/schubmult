from schubmult import *
from schubmult.symbolic import *
from sympy import pretty_print
import itertools

def hom_into_grass_ring(rc, grass_ring):
    if isinstance(rc, RCGraph):
        result = grass_ring.zero
        cem = RCGraph.full_CEM(rc.perm, len(rc))
        result += grass_ring.from_dict(cem.get(rc, {}))
        if rc.is_principal:
            for rcc, cem_dict in cem.items():
                if rcc.perm != rc.perm:
                    result += grass_ring.from_dict(cem_dict)
        assert expand(result.to_rc_graph_ring_element().polyvalue(Sx.genset) - rc.polyvalue(Sx.genset)) == S.Zero
        return result
    result = grass_ring.zero
    for rcc, coeff in rc.items():
        result += coeff * hom_into_grass_ring(rcc, grass_ring)
    return result

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)
    ring = GrassTensorAlgebra()
    r = RCGraphRing()

    ghost_rcs = {}
    for perm in perms:
        for rc, cem_dict in RCGraph.full_CEM(perm, n - 1).items():
            if rc.is_highest_weight and rc.perm != perm:
                ghost_rcs[(rc, perm)] = cem_dict
    ideal = {}
    for perm in perms:
        for rc, cem_dict in RCGraph.full_CEM(perm, n - 1).items():
            
            elem = ring.from_dict(cem_dict)
            for (rcc, perm2), cem_dict2 in ghost_rcs.items():
                elem2 = ring.from_dict(cem_dict2)
                try1 = (elem * elem2).to_rc_graph_ring_element().resize(n - 1)
                if try1 != r.zero:
                    assert len(try1) == 2, f"Expected exactly two terms in {try1}"
                    key, value = None, None
                    for rcc2, coeff in try1.items():
                        if coeff == -1:
                            key = rcc2.to_highest_weight()[0]
                        elif coeff == 1:
                            value = rcc2.to_highest_weight()[0]
                        else:
                            raise ValueError(f"Unexpected coefficient {coeff} in {try1}")
                    ideal[key] = ideal.get(key, set())
                    
                    ideal[key].add(value)
                    pretty_print(key)
                    print("Mapsto")
                    pretty_print(value)
                    # check equiv
                    keyn = key.to_highest_weight()[0]
                    valuen = value.to_highest_weight()[0]
                    keyn1, keyn2 = keyn.squash_decomp()
                    valuen1, valuen2 = valuen.squash_decomp()
                    assert keyn2.perm.inv != 0
                    assert valuen2.perm.inv != 0
                try2 = (elem2 * elem).to_rc_graph_ring_element().resize(n - 1)
                if try2 != r.zero:
                    assert len(try2) == 2, f"Expected exactly two terms in {try2}"
                    key, value = None, None
                    for rcc2, coeff in try2.items():
                        if coeff == -1:
                            key = rcc2.to_highest_weight()[0]
                        elif coeff == 1:
                            value = rcc2.to_highest_weight()[0]
                        else:
                            raise ValueError(f"Unexpected coefficient {coeff} in {try2}")
                    if key in ideal:
                        ideal[key].add(value)
                    else:
                        ideal[key] = {value}
                    pretty_print(key)
                    print("Mapsto")
                    pretty_print(value)
                    keyn = key.to_highest_weight()[0]
                    valuen = value.to_highest_weight()[0]
                    keyn1, keyn2 = keyn.squash_decomp()
                    valuen1, valuen2 = valuen.squash_decomp()
                    #assert keyn1 == valuen1 and keyn2 == valuen2, f"Failed equiv check for {key} and {value} with keyn {keyn} and valuen {valuen} and decomps {keyn1, keyn2} and {valuen1, valuen2}"
                    assert keyn2.perm.inv != 0
                    assert valuen2.perm.inv != 0
    print(f"Found {len(ideal)} elements in the ideal.")
                