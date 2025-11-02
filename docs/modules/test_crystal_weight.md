<!-- filepath: src/tests/test_crystal_weight.py -->

- `FunctionDef` — `_pad_left` (line 7)
- `FunctionDef` — `_add_weights` (line 12)
- `FunctionDef` — `rcgraph_crystal_weight_via_phi_eps` (line 18)
    - Recompute crystal weight using phi(i) - epsilon(i) for i = 1..(n-1)
- `FunctionDef` — `compute_weights_for_element` (line 42)
    - For mapping-like elements (with .items()) return set of weights of basis keys
- `FunctionDef` — `test_crystal_weight_matches_phi_eps_small_perms` (line 85)
- `FunctionDef` — `test_weights_for_ring_element` (line 107)
- `FunctionDef` — `test_weights_for_tensor_element` (line 123)
    - Construct a tensor ring element (tuple-key) and verify its weight equals