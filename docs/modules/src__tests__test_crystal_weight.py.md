# src/tests/test_crystal_weight.py



## _pad_left(t, length)



## _add_weights(w1, w2)



## rcgraph_crystal_weight_via_phi_eps(rc)

Recompute crystal weight using phi(i) - epsilon(i) for i = 1..(n-1)
and then cumulative sums as used in the codebase.

## compute_weights_for_element(elem)

For mapping-like elements (with .items()) return set of weights of basis keys
with nonzero coefficient. For RCGraph input, return single tuple.
For tuple keys (rc1, rc2, ...), compute elementwise sum of component weights.

## test_crystal_weight_matches_phi_eps_small_perms()



## test_weights_for_ring_element()



## test_weights_for_tensor_element()

Construct a tensor ring element (tuple-key) and verify its weight equals
the elementwise sum of component RCGraph weights.

