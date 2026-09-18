from schubmult import *

def cancel_forest(elem):
    the_key_dict = {}
    new_elem = 0
    for key, value in elem.items():
        if key.forest_invariant in the_key_dict:
            new_elem += value * elem.ring(the_key_dict[key.forest_invariant])
        else:
            new_elem += value * elem.ring(key)
            the_key_dict[key.forest_invariant] = key
    return new_elem

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    r = RCGraphRing()
    for perm in Permutation.all_permutations(n):
        for length in range(max(2,len(perm.trimcode)), n):
            #for rc in RCGraph.all_hw_rcs(perm, length):
            #    rc_cut, row = rc.vertical_cut(length - 1)
            ring_elem = r.from_free_algebra_element(ASx(perm, length))
            ring_elem = cancel_forest(ring_elem)
            assert ring_elem.to_free_algebra_element() == ASx(perm, length), f"Error: after canceling forest invariant duplicates in {perm.trimcode}, got {ring_elem.to_free_algebra_element()} instead of {ASx(perm, length)}"
            
            assert all(v >= 0 for _, v in ring_elem.items()), f"Negative coefficient in coproduct of {perm.trimcode}, got {ring_elem}"
