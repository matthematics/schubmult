from schubmult import *
from schubmult.rings.polynomial_algebra import *
from schubmult.symbolic.common_polys import lascoux_poly

def _bwify_key(key):
    new_key = [next(iter(RCGraph.elem_sym_rcs(len(k.perm_word), k.perm.max_descent, weight=k.length_vector))) for k in key]
    return br.make_key(new_key, key.size)

if __name__ == "__main__":
    import sys
    import itertools
    from sympy import pretty_print
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    br = BoundedRCFactorAlgebra()
    bw = BoundedWCFactorAlgebra()
    dct = {}
    rc_ring = RCGraphRing()
    wc_ring = WCGraphRing()
    perms_dd = Permutation.all_permutations(n * 2 - 1)
    def _morphism(rc):
        cheat_elem = ASx(rc.perm, len(rc)).change_basis(GrothendieckBasis)
        ret = wc_ring.zero
        for (perm, length), coeff in cheat_elem.items():
            if coeff == 0:
                continue
            #for wc in WCGraph.all_wc_graphs(perm, length, weight=rc.length_vector):
            g = bw.full_groth_elem(perm, len(perm), 1)
            for cheat_key, coeff in g.items():
                rcc = br.key_to_rc_graph(_bwify_key(cheat_key)).resize(len(rc))
                if rcc == rc:
                    wc = bw.key_to_wc_graph(cheat_key).resize(len(rc))
                    ret += coeff * wc_ring(wc)
            # key = next(iter(br.from_rc_graph))
            # if br.key_to_rc_graph(_bwify_key(br.make_key([wc], 1))).resize(len(rc)).almosteq(rc):
            #     ret += coeff * wc_ring(wc)
        return ret

    def morphism(rc_ring_elem):
        ret = 0
        for rc, coeff in rc_ring_elem.items():
            ret += coeff * _morphism(rc)
        return ret

    # for perm in perms_dd:
    #     for rc in RCGraph.all_rc_graphs(perm, 2 *n - 1):
    #         print(f"perm={perm.trimcode}")
            
            
    #         #print(f"    rc={rc}  wc_graphs={dct[rc]}")
    comps = [perm.pad_code(n - 1) for perm in perms]
    
    for comp in comps:
        if sum(comp) == 0:
            continue
        key_result = 0
        actual_lascoux = KeyPoly.from_expr(lascoux_poly(comp, Sx.genset), length=n-1)
        principal = WCGraph(RCGraph.principal_rc(uncode(comp), n - 1))
        # for wc in WCGraph.all_wc_graphs(uncode(comp), n - 1):
        #     if wc.hecke_invariant[0] == principal.hecke_invariant[0]:
        #         key_result += wc.polyvalue(Sx.genset, beta=1)
        
        g = bw.full_groth_elem(uncode(comp), n + 2, 1)
        # #for wc in {wcc for wcc in WCGraph.all_wc_graphs(uncode(comp), n - 1) if wcc.hecke_invariant[0] == WCGraph(RCGraph.principal_rc(uncode(comp), n - 1)).hecke_invariant[0]}:
        # principal = WCGraph(RCGraph.principal_rc(uncode(comp), n - 1))
        for key, coeff in g.items():
            wc = bw.key_to_wc_graph(key).resize(n - 1)
            if wc.perm == uncode(comp) and wc.strong_hecke_invariant == principal.strong_hecke_invariant:
                rc = br.key_to_rc_graph(_bwify_key(key)).resize(n - 1)
                if rc.is_highest_weight:
                    key_result += KeyPoly(*rc.extremal_weight)
                    #key_result += rc.polyvalue(Sx.genset, crystal=True)
                #key_result += wc.polyvalue(Sx.genset, beta=1)
                #key_result += wc.length_vector.count(1) * morphism(rc_ring(rc_ring.from_wc_graph(wc)))
        #key_result = KeyPoly.from_expr(key_result, length=n-1)
        assert key_result.almosteq(actual_lascoux), f"Mismatch for {comp}: {key_result} != {actual_lascoux}\n{key_result-actual_lascoux}"
            