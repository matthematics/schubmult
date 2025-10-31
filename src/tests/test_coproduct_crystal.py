import pytest

from schubmult import ASx, Permutation
from schubmult.schub_lib.rc_graph import RCGraph
from schubmult.schub_lib.rc_graph_ring import RCGraphRing


def _normalize_fa_coprod(fa_coprod):
    """
    Normalize the various possible key shapes from ASx(...).coproduct()
    into a mapping (perm_u, perm_v) -> multiplicity.
    """
    out = {}
    for key, mult in fa_coprod.items():
        # Try to unpack (perm_u, perm_v) directly
        try:
            pu, pv = key
        except Exception:
            # Unexpected shape: skip it
            continue

        # If items are (perm, length) pairs, unwrap
        if isinstance(pu, tuple) and len(pu) >= 1 and hasattr(pu[0], "trimcode"):
            pu = pu[0]
        if isinstance(pv, tuple) and len(pv) >= 1 and hasattr(pv[0], "trimcode"):
            pv = pv[0]

        # If items are RCGraph-like, use .perm
        if hasattr(pu, "perm"):
            pu = pu.perm
        if hasattr(pv, "perm"):
            pv = pv.perm

        out[(pu, pv)] = out.get((pu, pv), 0) + mult
    return out


def _project_tensor_to_perms(tensor_elem):
    """
    Project a tensor-module element (keys are (RCGraph, RCGraph)) to
    (perm_u, perm_v) -> multiplicity by summing coefficients of pairs
    that share permutations.
    """
    out = {}
    for (L, R), mult in tensor_elem.items():
        pu = getattr(L, "perm", L)
        pv = getattr(R, "perm", R)
        out[(pu, pv)] = out.get((pu, pv), 0) + mult
    return out


@pytest.mark.parametrize("n", [3, 4])
def test_coproduct_on_basis_principal_matches_ASx(n):
    """
    For principal RC graphs (constructed with RCGraph.principal_rc), the
    RCGraphRing.coproduct_on_basis projection to permutation pairs should
    match the permutation-level ASx(...).coproduct() multiplicities.
    This validates the LR/crystal matching path for small n.
    """
    length = n - 1
    ring = RCGraphRing()

    perms = Permutation.all_permutations(n)
    for perm in perms:
        rc = RCGraph.principal_rc(perm, length)

        # permutation-level Schubert coproduct
        fa = ASx(perm, length).coproduct()
        fa_norm = _normalize_fa_coprod(fa)

        # RCGraph-level coproduct (tensor of RCGraph x RCGraph)
        tr = ring.coproduct_on_basis(rc)
        tr_proj = _project_tensor_to_perms(tr)

        # Compare multiplicities for union of keys
        all_keys = set(fa_norm.keys()) | set(tr_proj.keys())
        for k in all_keys:
            a = fa_norm.get(k, 0)
            b = tr_proj.get(k, 0)
            assert a == b, (
                f"Coproduct mismatch for perm {perm} (n={n}) on pair {k}:\n"
                f"  ASx multiplicity = {a}\n"
                f"  RCGraph coproduct  = {b}\n"
            )