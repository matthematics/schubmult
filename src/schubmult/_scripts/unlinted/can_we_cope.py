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
    from schubmult.rings.rc_graph_ring import RCGraphRing

    res = S.Zero
    for rc, coeff in rc_elem.items():
        print(F"{rc=}")
        print(f"{coeff=}")
        cprd = rc_coproduct(rc)
        res += coeff * cprd
    print(f"{rc_elem=}")
    print(f"{res=}")
    return res

if __name__ == "__main__":
    import sys
    import itertools

    n = int(sys.argv[1])

    # map RC to Schub and weight, coproduct, recombine

    perms = Permutation.all_permutations(n)
    
    for perm1, perm2 in itertools.product(perms, repeat=2):
        length = max(len(perm1.trimcode), len(perm2.trimcode))
        for rc1,rc2 in itertools.product(RCGraph.all_rc_graphs(perm1, length), RCGraph.all_rc_graphs(perm2, length)):
            rc_elem = R(rc1) * R(rc2)
            cprd1 = rc_coproduct(rc1)
            cprd2 = rc_coproduct(rc2)
            cprd_total = cprd1 * cprd2

            rc_cprd_total = rc_ring_coproduct(rc_elem)
            try:
                assert all(v == S.Zero for v in (rc_cprd_total - cprd_total).values()), f"Error: RC graph ring element coproduct mismatch for permutations {perm1} and {perm2}, {rc1}, {rc2}\nComputed:\n{rc_cprd_total}\nExpected:\n{cprd_total}"
            except AssertionError as e:
                print(f"Computed:\n{rc_cprd_total}\nExpected:\n{cprd_total}")
                print(rc_cprd_total - cprd_total)
                sys.exit(1)
            
                
                