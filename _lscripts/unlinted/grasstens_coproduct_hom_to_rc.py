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

    for perm1, perm2 in itertools.product(perms, repeat=2):
        prd = Sx(perm1) * Sx(perm2)
        rprd = r.zero
        for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1, n - 1), RCGraph.all_rc_graphs(perm2, n - 1)):
            # elem1 = hom_into_grass_ring(r(rc1)*r(rc2), ring).to_rc_graph_ring_element()
            # elem2 = hom_into_grass_ring(rc1, ring).dual_product(hom_into_grass_ring(rc2, ring)).to_rc_graph_ring_element()
            # assert elem1.almosteq(elem2), f"Failed for {rc1} and {rc2} with product {r(rc1) * r(rc2)} {elem1=} {elem2=}"
            rprd += (hom_into_grass_ring(rc1, ring) * hom_into_grass_ring(rc2, ring)).to_rc_graph_ring_element()
        for rcc, coeff in rprd.items():
            assert prd.get(rcc.perm, 0) == coeff, f"Failed for {rc1} and {rc2} with product {prd} {rprd=} {prd=} {rcc.perm=}"
            #cprd = hom_into_grass_ring(r(rc1) * r(rc2), ring).coproduct()
            # grass_elem2 = hom_into_grass_ring(r(rc1) * r(rc2), ring)
            # for key, coeff in grass_elem2.items():
            #     assert coeff == 1
            #     cprd = ring(key).coproduct()
            #     pretty_print(cprd)
            #     for (key1, key2), coeff2 in cprd.items():
            #         rc1_elem = next(iter(ring(key1).to_rc_graph_ring_element()))
            #         rc2_elem = next(iter(ring(key2).to_rc_graph_ring_element()))
            #         assert coeff2 == 1
            #         if len(rc1_elem) == len(rc1) and len(rc2_elem) == len(rc2):
            #             assert rc1_elem == rc1, f"Failed for {key1} with coproduct {cprd} {rc1_elem=} {rc1=}"
            #             assert rc2_elem == rc2, f"Failed for {key2} with coproduct {cprd} {rc2_elem=} {rc2=}"
                    # if all(len(rcc) == len(rc1) for rcc in rc1_elem.keys()) and all(len(rcc) == len(rc2) for rcc in rc2_elem.keys()):
                    #     assert rc1_elem.almosteq(r(rc1)) and rc2_elem.almosteq(r(rc2)), f"Failed for {key1} and {key2} with coproduct {cprd} {rc1_elem=} {rc2_elem=} {r(rc1)=} {r(rc2)=}"
                    #assert rc2_elem.almosteq(r(rc2)), f"Failed for {key2} with coproduct {cprd}"
    #     # #for rc in RCGraph.all_rc_graphs(perm):
    #     # cem = RCGraph.full_CEM(perm, n - 1)
    #     # for rc, cem_dict in cem.items():
    #     #     elem = ring.from_dict(cem_dict)
    #     #     cprd = elem.coproduct()
    #     #     print("RC:")
    #     #     pretty_print(rc)
    #     #     print("cprd:")
    #     #     pretty_print(cprd)
    # Collect kernel generators: differences of tensor keys that map to the same RC graph
    # representatives = {}
    # kernel_elems = []
    # for perm in perms:
    #     cem = RCGraph.full_CEM(perm, n - 1, tuple(range(n, 0, -1)))
    #     for rc, cem_dict in cem.items():
    #         if not rc.is_highest_weight:
    #             continue
    #         elem = ring.from_dict(cem_dict)
    #         if rc.perm != perm:
    #             kernel_elems.append(elem)
    #         else:
    #             representatives[rc] = representatives.get(rc, [])
    #             representatives[rc].append(elem)

    # # The real kernel elements come from linear dependencies:
    # # find pairs of GrassTensorAlgebra elements with same image under pi
    # # for key1, key2 in pairs_with_same_rc_image:
    # #     k = ring(key1) - ring(key2)  # kernel element
    # #     for basis_key in all_basis_keys:
    # #         assert (ring(basis_key) * k).to_rc_graph_ring_element().almosteq(r.zero)
    # #         assert (k * ring(basis_key)).to_rc_graph_ring_element().almosteq(r.zero)
    # for k in kernel_elems:
    #     for rc, rep_list in representatives.items():
    #         for rep in rep_list:
    #             assert (rep * k).to_rc_graph_ring_element().almosteq(r.zero), f"Failed right kernel test for {k} and representative {(rep * k).to_rc_graph_ring_element()} of {rc}"
    #             assert (k * rep).to_rc_graph_ring_element().almosteq(r.zero), f"Failed left kernel test for {k} and representative {(k * rep).to_rc_graph_ring_element()} of {rc}"