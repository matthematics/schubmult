# src/scripts/nr_verification.py



## _normalize_fa_coprod(fa_coprod)

Normalize the various possible key shapes from ASx(...).coproduct()
into a mapping (perm_u, perm_v) -> multiplicity.

## _project_tensor_to_perms(tensor_elem)

Project a tensor-module element (keys are (RCGraph, RCGraph)) to
(perm_u, perm_v) -> multiplicity by summing coefficients of pairs
that share permutations.

## main()



