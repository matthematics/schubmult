# src/tests/test_coproduct_crystal.py



## _normalize_fa_coprod(fa_coprod)

Normalize the various possible key shapes from ASx(...).coproduct()
into a mapping (perm_u, perm_v) -> multiplicity.

## _project_tensor_to_perms(tensor_elem)

Project a tensor-module element (keys are (RCGraph, RCGraph)) to
(perm_u, perm_v) -> multiplicity by summing coefficients of pairs
that share permutations.

## test_coproduct_on_basis_principal_matches_ASx(n)

For principal RC graphs (constructed with RCGraph.principal_rc), the
RCGraphRing.coproduct_on_basis projection to permutation pairs should
match the permutation-level ASx(...).coproduct() multiplicities.
This validates the LR/crystal matching path for small n.

