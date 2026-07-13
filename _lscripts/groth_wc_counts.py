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
    
    for u in perms:
        g = bw.full_groth_elem(u, n)
        dct = {}
        for key, coeff in g.items():
            wc = bw.key_to_wc_graph(key)
            if wc.perm == u:
                dct[wc] = dct.get(wc, set())
                dct[wc].add(key)
        
        for wc, keys in dct.items():
            assert len(keys) == 1, f"Multiple keys for {wc}: {keys}"
    