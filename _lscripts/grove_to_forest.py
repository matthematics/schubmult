from schubmult import *

bw = BoundedWCFactorAlgebra()
br = BoundedRCFactorAlgebra()

def _bwify_key(key):
    new_key = [next(iter(RCGraph.elem_sym_rcs(len(k.perm_word), k.perm.max_descent, weight=k.length_vector))) for k in key]
    return br.make_key(new_key, key.size)

def hom_bw_to_br(bw_elem):
    """Map a BoundedWCFactorAlgebra element to a BoundedRCFactorAlgebra element."""
    ret = 0
    for key, coeff in bw_elem.items():
        new_key = _bwify_key(key)
        ret += coeff * br(new_key)
    return ret

def _g_grove_extractor(rc_elem, indexes):
    """The operator G_i on rc graphs, which is the composition of the quasi-shift and trim-descent operations."""
    result = f.zero
    if len(indexes) > 0:
        desc, is_left_child = indexes[-1]      
        for rc, coeff in rc_elem.items():
            forest = weak_composition_to_indfor(rc.forest_weight)
            if desc not in forest.trim_descents:
                continue
            if is_left_child:
                result += coeff * f(rc).forest_trim(desc)
            else:
                result += coeff * (f(rc).forest_trim(desc) - f(rc).quasi_shift(desc) + f(rc).quasi_shift(desc + 1))
        if len(indexes) > 1:
            return _g_grove_extractor(result, indexes=indexes[:-1])
    else:
        result = rc_elem
    result_values = [v for k, v in result.items() if k.perm.inv == 0]
    return sum(result_values)

if __name__ == "__main__":
    import sys
    import itertools
    from schubmult.rings.polynomial_algebra import *
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    
    for u in perms:
        g = bw.full_groth_elem(u, n)
        dct = {}
        grovy = GrovePoly(*u.pad_code(n - 1))
        forestify = 0
        for key, coeff in g.items():
            wc = bw.key_to_wc_graph(key).resize(n - 1)
            
            if wc.grove_weight == u.pad_code(n - 1):
                assert coeff > 0, f"Unexpected negative coefficient {coeff} for {key=}, {wc=}"
                rc = br.key_to_rc_graph(_bwify_key(key)).resize(n - 1)
                if rc.forest_weight == rc.length_vector:
                    forestify += coeff * ForestPoly(*rc.forest_weight)
        assert grovy.change_basis(ForestPolyBasis).almosteq(forestify), f"Mismatch for {u}: {grovy.change_basis(ForestPolyBasis)} != {forestify}"
        print(f"Grab the bag {u}")
        # for wc, keys in dct.items():
        #     assert len(keys) == 1, f"Multiple keys for {wc}: {keys}"
    