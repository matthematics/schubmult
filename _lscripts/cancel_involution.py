"""Probe the sign-cancellation structure of the signed key-crystal model.

For each permutation we build the highest-weight crystal dictionary exactly as
``groth_key_decomp.main`` does, then group the crystals by their key label.
Within each key-label group we print every distinct ``CrystalGraphTensor``
(spatula) together with:

  * its signed coefficient,
  * its crystal weight,
  * a factorization signature: for each elem-sym RC factor, the pair
    ``(perm_word, length_vector, excess)``.

The goal is to see what distinguishes the two members of a cancelling +/- pair
so that a canonical sign-reversing involution can be defined on the factor data.
"""

import sys
from collections import defaultdict
from math import prod

from schubmult import CrystalGraphTensor, Gx, Permutation, Sx
from schubmult.symbolic import expand, sympify
from schubmult.rings.polynomial_algebra import KeyPolyBasis, PolynomialAlgebra

sys.argv = sys.argv  # keep argv
import importlib.util

_spec = importlib.util.spec_from_file_location(
    "gkd", __file__.replace("cancel_involution.py", "groth_key_decomp.py")
)
gkd = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(gkd)

KeyPoly = PolynomialAlgebra(KeyPolyBasis(Gx.genset))
beta = Gx._beta


def factor_signature(spatula):
    """Return a tuple describing each elem-sym RC factor of the tensor."""
    sig = []
    for factor in spatula:
        pw = getattr(factor, "perm_word", None)
        lv = getattr(factor, "length_vector", None)
        exc = getattr(factor, "excess", None)
        sig.append((tuple(pw) if pw is not None else None,
                    tuple(lv) if lv is not None else None,
                    exc))
    return tuple(sig)


def key_label(spatula):
    spanko = sum(prod([e.polyvalue(Sx.genset) for e in kk]) for kk in spatula.full_crystal)
    kp = KeyPoly.from_expr(spanko)
    labels = [p for p, v in kp.items() if sympify(v).expand() != 0]
    return tuple(labels[0]) if labels else None


def build_hw_dict(perm, length):
    tagged = gkd.untagged_groth_elem(perm, length)
    hw_dict = {}
    for key, coeff in tagged.items():
        if key.is_highest_weight:
            the_excess = sum(k.excess for k in key)
            spatula = CrystalGraphTensor(*[k._convert_elem_rc() for k in key])
            hw_dict[spatula] = hw_dict.get(spatula, 0) + beta ** the_excess * coeff
    return hw_dict


def is_principal_wide(word):
    """A factor is a principal wide elem sym if its perm_word is (1,2,...,k), k>=2."""
    return word is not None and len(word) >= 2 and tuple(word) == tuple(range(1, len(word) + 1))


def num_principal_wide(multiset):
    return sum(1 for f in multiset if is_principal_wide(f[0]))


def main(n):
    length = n + 2
    all_ok = True
    for perm in Permutation.all_permutations(n):
        if perm.inv == 0:
            continue
        hw_dict = build_hw_dict(perm, length)
        # aggregate by (key label, UNORDERED factor multiset) to remove tensor-order artifacts
        agg = defaultdict(lambda: 0)
        cwmap = {}
        for spatula, coeff in hw_dict.items():
            kl = key_label(spatula)
            sig = factor_signature(spatula)
            multiset = tuple(sorted(sig))
            akey = (kl, multiset)
            agg[akey] = agg[akey] + sympify(coeff)
            cwmap[akey] = spatula.crystal_weight

        # CANONICAL POSITIVE MODEL: keep only factorization types with NO principal
        # wide factor. Sum their (positive) contributions as key polynomials.
        survivor_poly = 0
        survivor_terms = []
        for (kl, multiset), coeff in agg.items():
            c = expand(coeff)
            if c == 0:
                continue
            if num_principal_wide(multiset) == 0:
                survivor_poly += c * KeyPoly(kl).expand()
                survivor_terms.append((kl, multiset, c))

        target = Gx(perm).expand()
        ok = expand(survivor_poly - target) == 0
        # are all survivor coeffs positive at beta=1?
        pos = all(sympify(c).subs(beta, 1) > 0 for _, _, c in survivor_terms)
        all_ok = all_ok and ok and pos
        flag = "OK " if (ok and pos) else "FAIL"
        print(f"[{flag}] perm={perm.trimcode}  survivors={len(survivor_terms)}  reproduces_Gx={ok}  all_positive={pos}")
        if not ok:
            print(f"      residual = {expand(survivor_poly - target)}")
    print(f"\nALL_OK={all_ok}")


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    main(n)
