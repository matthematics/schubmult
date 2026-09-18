from schubmult import RCGraph
from schubmult import RCGraphRing
from sympy import pretty_print

from schubmult import ASx, Permutation


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
# LR rule verification script

import os
import shutil
import sys


def main():

    n = int(sys.argv[1])

    length = n - 1
    ring = RCGraphRing()

    perms = Permutation.all_permutations(n)
    suc_count = 0
    tot = 0
    for perm in perms:
        for rc in RCGraph.all_rc_graphs(perm, n-1):

            # permutation-level Schubert coproduct
            fa = ASx(perm, length).coproduct()
            fa_norm = _normalize_fa_coprod(fa)

            # RCGraph-level coproduct (tensor of RCGraph x RCGraph)
            tr = ring.coproduct_on_basis(rc)
            if tr is None:
                print(f"Skipping perm {perm}: coproduct_on_basis returned None")
                continue
            tr_proj = _project_tensor_to_perms(tr)


            # Compare multiplicities for union of keys
            all_keys = set(fa_norm.keys()) | set(tr_proj.keys())
            success = True
            for k in all_keys:
                a = fa_norm.get(k, 0)
                b = tr_proj.get(k, 0)
                try:
                    assert a == b, (
                        f"Coproduct mismatch for perm {perm} (n={n}) on pair {k}:\n"
                        f"  ASx multiplicity = {a}\n"
                        f"  RCGraph coproduct  = {b}\n"
                    )
                except AssertionError as e:
                    print(e)
                    success = False
                    break
            pretty_print(tr)
            if success:
                print(f"Coproduct match for perm {perm} (n={n})")
                suc_count += 1
            tot += 1
    print(f"Successful coproduct matches: {suc_count} / {tot} for n={n}")

if __name__ == "__main__":
    main()
