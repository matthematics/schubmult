from schubmult import *
from functools import cache


def _bwify_key(key):
    new_key = [next(iter(RCGraph.elem_sym_rcs(len(k.perm_word), k.perm.max_descent, weight=k.length_vector))) for k in key]
    return br.make_key(new_key, key.size)

@cache
def _get_groth_elem(perm, length):
    return bw.full_groth_elem(perm, length, 1)

@cache
def _cheat_elem(perm, length):
    return ASx(perm, length).change_basis(GrothendieckBasis)

_cheat_stash = set()

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
    
    @cache
    def _morphism(rc):
        global _cheat_stash
        length = len(rc)
        cheat_elem = _cheat_elem(rc.perm, length)
        ret = wc_ring.zero
        
        
        for (perm, length), coeff1 in cheat_elem.items():
            if coeff1 == 0:
                continue
            g = _get_groth_elem(perm, length)
            for cheat_key, coeff2 in g.items():
                if coeff2 == 0:
                    continue
                rcc = br.key_to_rc_graph(_bwify_key(cheat_key)).resize(len(rc))
                if rcc == rc:
                    wc = bw.key_to_wc_graph(cheat_key).resize(len(rc))
                    ret += coeff2 * wc_ring(wc)
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
    for perm1, perm2 in itertools.combinations(perms, 2):
        if perm1.inv == 0 or perm2.inv == 0:
            continue
        for rc1 in RCGraph.all_rc_graphs(perm1):
            m1 = morphism(rc_ring(rc1))
            for rc2 in RCGraph.all_rc_graphs(perm2):
                wc_prod = m1 * morphism(rc_ring(rc2))
                wc_rc_prod = morphism(rc_ring(rc1) * rc_ring(rc2))
                try:
                    assert wc_rc_prod.almosteq(wc_prod), f"Mismatch for {rc1}, {rc2}: {wc_prod-wc_rc_prod}"
                except AssertionError as e:
                    pretty_print(wc_prod)
                    pretty_print(wc_rc_prod)
                    raise e
                print("ASIFN")