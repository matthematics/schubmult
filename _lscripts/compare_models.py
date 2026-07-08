"""Confirm the signed CEM key model and the positive RC-crystal key model agree.

- signed model: ``groth_key_decomp.untagged_groth_elem`` -> collapse highest
  weight crystals to key labels (signs present, cancel by linear independence).
- positive model: group ``RCGraph.all_rc_graphs`` of each Schubert term of
  ``WCGraph.groth_to_schub(perm, beta)`` into Demazure crystals by highest
  weight; each fiber is one ``Key(extremal_weight)`` with a positive beta-power.

We print the key expansion from each side and assert they match.
"""

import importlib.util
import sys
from collections import defaultdict
from math import prod

from schubmult import CrystalGraphTensor, Gx, Permutation, RCGraph, Sx, WCGraph
from schubmult.symbolic import expand, sympify
from schubmult.rings.polynomial_algebra import KeyPolyBasis, PolynomialAlgebra

_spec = importlib.util.spec_from_file_location(
    "gkd", __file__.replace("compare_models.py", "groth_key_decomp.py")
)
gkd = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(gkd)

x = Sx.genset
KeyPoly = PolynomialAlgebra(KeyPolyBasis(x))
beta = Gx._beta


def _strip(label):
    if label is None:
        return None
    t = tuple(label)
    while t and t[-1] == 0:
        t = t[:-1]
    return t


def signed_key_expansion(perm, length):
    """Collapse the signed hw-crystal model to {key_label: coeff}."""
    tagged = gkd.untagged_groth_elem(perm, length)
    out = defaultdict(lambda: 0)
    for key, coeff in tagged.items():
        if not key.is_highest_weight:
            continue
        the_excess = sum(k.excess for k in key)
        spatula = CrystalGraphTensor(*[k._convert_elem_rc() for k in key])
        spanko = sum(prod([e.polyvalue(x) for e in kk]) for kk in spatula.full_crystal)
        kp = KeyPoly.from_expr(spanko)
        labels = [p for p, v in kp.items() if sympify(v).expand() != 0]
        kl = _strip(labels[0]) if labels else None
        out[kl] = out[kl] + beta ** the_excess * sympify(coeff)
    return {k: expand(v) for k, v in out.items() if expand(v) != 0}


def rc_key_expansion(perm, length):
    """Positive RC-crystal key model: {key_label: coeff}."""
    out = defaultdict(lambda: 0)
    for perm2, coeff in WCGraph.groth_to_schub(perm, beta).items():
        by_hw = defaultdict(list)
        for rc in RCGraph.all_rc_graphs(perm2, length):
            by_hw[rc.to_highest_weight()[0]].append(rc)
        for hw, _ in by_hw.items():
            kl = _strip(hw.extremal_weight)
            out[kl] = out[kl] + sympify(coeff)
    return {k: expand(v) for k, v in out.items() if expand(v) != 0}


def main(n):
    length = n + 2
    all_ok = True
    for perm in Permutation.all_permutations(n):
        if perm.inv == 0:
            continue
        signed = signed_key_expansion(perm, length)
        rc = rc_key_expansion(perm, length)
        keys = set(signed) | set(rc)
        match = all(expand(signed.get(k, 0) - rc.get(k, 0)) == 0 for k in keys)
        rc_positive = all(sympify(v).subs(beta, 1) > 0 for v in rc.values())
        all_ok = all_ok and match and rc_positive
        flag = "OK " if (match and rc_positive) else "FAIL"
        print(f"[{flag}] perm={perm.trimcode}  keys={len(rc)}  models_agree={match}  rc_positive={rc_positive}")
        if not match:
            for k in sorted(keys, key=str):
                s, r = expand(signed.get(k, 0)), expand(rc.get(k, 0))
                if s != r:
                    print(f"      key={k}: signed={s}  rc={r}")
    print(f"\nALL_OK={all_ok}")


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    main(n)
