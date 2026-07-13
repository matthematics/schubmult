from schubmult import *

bw = BoundedWCFactorAlgebra()
br = BoundedRCFactorAlgebra()

def hom_bw_to_br(bw_elem):
    """Map a BoundedWCFactorAlgebra element to a BoundedRCFactorAlgebra element."""
    ret = 0
    for key, coeff in bw_elem.items():
        new_key = [next(iter(RCGraph.elem_sym_rcs(len(k.perm_word), k.perm.max_descent, weight=k.length_vector))) for k in key]
        ret += coeff * br(br.make_key(new_key, key.size))
    return ret


if __name__ == "__main__":
    import sys
    import itertools
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    for u, v in itertools.product(perms, repeat=2):
        g1 = bw.full_groth_elem(u, n)
        g2 = bw.full_groth_elem(v, n)
        prod = g1 * g2
        prd1 = hom_bw_to_br(prod).to_rc_graph_ring_element()
        prd2 = (hom_bw_to_br(g1) * hom_bw_to_br(g2)).to_rc_graph_ring_element()
        if not prd1.almosteq(prd2):
            print(f"Mismatch for {u}, {v}: {prd1} != {prd2}")
    