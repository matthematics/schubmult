from schubmult import *
from schubmult.symbolic import S, expand

FA = FreeAlgebra(WordBasis)
T = ASx @ FA

def rc_coproduct(rc: RCGraph) -> dict:
    """
    Compute the coproduct of an RC graph into pairs of RC graphs.

    Args:
        rc (RCGraph): The RC graph to coproduct.

    Returns:
        dict: A dictionary mapping ((perm1, len1), (perm2, len2)) to coefficient.
    """
    schub_elem = ASx(rc.perm, len(rc))
    weight_elem = FA(*rc.length_vector)
    cprd1 = schub_elem.coproduct()
    cprd2 = weight_elem.coproduct()
    total_coprod = (T@T).zero
    length = len(rc)
    for ((perm1, _), (perm2, _)), coeff1 in cprd1.items():
        for (w1, w2), coeff2 in cprd2.items():
            total_coprod += coeff1 * ASx(perm1, length).change_basis(WordBasis).get(w1, S.Zero) * ASx(perm2, length).change_basis(WordBasis).get(w2, S.Zero)*(T@T).ext_multiply(T.ext_multiply(ASx(perm1, length), FA(*w1)), T.ext_multiply(ASx(perm2, length), FA(*w2)))
    return total_coprod

R = RCGraphRing()

def rc_ring_coproduct(rc_elem: RCGraphRingElement) -> RCGraphRingElement:
    """
    Compute the coproduct of an RC graph ring element.

    Args:
        rc_elem (RCGraphRingElement): The RC graph ring element to coproduct.

    Returns:
        RCGraphRingElement: The coproduct as an RC graph ring element.
    """
    from schubmult.rings.combinatorial.rc_graph_ring import RCGraphRing

    res = S.Zero
    for rc, coeff in rc_elem.items():
        print(F"{rc=}")
        print(f"{coeff=}")
        cprd = rc_coproduct(rc)
        res += coeff * cprd
    print(f"{rc_elem=}")
    print(f"{res=}")
    return res

class RCDQRing(TensorRing):
    def __init__(self):
 
        super().__init__(RCGraphRing(), ASx)

    def mul(self, a, b):
        result = super().mul(a, b)
        result = super().from_dict({(k1, k2): v for (k1, k2), v in result.items() if v != 0 and k1.perm == k2[0] and len(k1) == k2[1]})
        return result

    def from_dict(self, dct):
        elem = {(k1, k2): v for (k1, k2), v in dct.items() if k1.perm == k2[0] and len(k1) == k2[1]}
        return super().from_dict(elem)

    def __call__(self, key_rc):
        return super().__call__((key_rc, (key_rc.perm, len(key_rc))))

def fa_to_rc_ring(w, w2, n):
    fa_elem = FA(*w)
    schub_elem = fa_elem.change_basis(SchubertBasis)
    rc_ring = RCGraphRing()
    ret = rc_ring.zero
    for (perm, length), coeff in schub_elem.items():
        if len(perm) > n:
            continue
        for rc in RCGraph.all_rc_graphs(perm, length, weight=w2):
            ret += rc.polyvalue(Sx.genset) * rc_ring(rc)
    return ret

if __name__ == "__main__":
    import sys
    import itertools

    n = int(sys.argv[1])

    # map RC to Schub and weight, coproduct, recombine

    perms = Permutation.all_permutations(n)
    perms_small = Permutation.all_permutations(n - 1)
    # for perm1, perm2 in itertools.product(perms, repeat=2):
    #     length = max(len(perm1.trimcode), len(perm2.trimcode))
    #     for rc1,rc2 in itertools.product(RCGraph.all_rc_graphs(perm1, length), RCGraph.all_rc_graphs(perm2, length)):
    #         # are distinct
    ring = RCDQRing()
    
    for perm in perms:
        
        A = ASx(perm, len(perm.trimcode)).change_basis(WordBasis)
        # res = (ring@ring).zero
        ring = RCGraphRing()
        T = ring @ ring
        test_elem = T.zero
        
        for w, coeff in A.items():
            #res += coeff * word_orc(*w)
            cprd = FA(*w).coproduct()
            code_prod = FA(*perm.trimcode).coproduct()
            for w1, w2 in cprd.keys():
                for c1, c2 in code_prod.keys():
                    test_elem += coeff * T.ext_multiply(fa_to_rc_ring(w1, c1, 100), fa_to_rc_ring(w2, c2, 100))
        #res = (ring@ring).from_dict({k: v for k, v in res.items() if v != 0 and len(k[0][1][0])<=n and len(k[1][1][0])<=n})
        print("------------------------------{perm}----------------------------")
        print("TEST:")
        print(test_elem)
        
            
                
                