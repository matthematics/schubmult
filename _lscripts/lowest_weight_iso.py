"""Anchor RC/WC crystals on the LOWEST (extremal) weight, per user's definition.

User's definition of the identifying (extremal) weight alpha of a crystal:
  1. alpha sorts (descending) to the highest weight lambda, and
  2. alpha is lex-minimal among all weights of the crystal with that property.

This alpha is the composition indexing the key kappa_alpha.

Step 1 (this script): verify the library's ``RCGraph.extremal_weight`` (a cumsum
dominance-min) equals the user's lex-min-weight alpha, computed directly from the
crystal's weight multiset.

Step 2: build the DIRECT anchor isomorphism RC-crystal <-> WC-crystal keyed on
alpha, and report the anchor elements (RC lowest-weight-demazure vs the WC element
carrying weight alpha).
"""

import sys
from collections import defaultdict

from schubmult import Gx, Permutation, RCGraph, Sx, WCGraph

beta = Gx._beta


def descending(t):
    return tuple(sorted(t, reverse=True))


def user_alpha(weights, lam):
    """Lex-min weight whose descending sort equals lam."""
    cands = [w for w in weights if descending(w) == tuple(lam)]
    return min(cands) if cands else None


def rc_crystal_data(perm, length):
    """Yield (perm2, hw_rc, crystal_list, lib_alpha, weights) for each RC Demazure crystal."""
    for perm2, coeff in WCGraph.groth_to_schub(perm, beta).items():
        by_hw = defaultdict(list)
        for rc in RCGraph.all_rc_graphs(perm2, length):
            by_hw[rc.to_highest_weight()[0]].append(rc)
        for hw, crystal in by_hw.items():
            full = list(hw.full_crystal)
            weights = [tuple(c.length_vector) for c in full]
            lam = tuple(hw.length_vector)
            yield perm2, coeff, hw, full, tuple(hw.extremal_weight), weights, lam


def main(n):
    length = n + 2
    mismatches = 0
    total = 0
    for perm in Permutation.all_permutations(n):
        if perm.inv == 0:
            continue
        for perm2, coeff, hw, full, lib_alpha, weights, lam in rc_crystal_data(perm, length):
            total += 1
            ua = user_alpha(weights, lam)
            # strip trailing zeros for comparison
            def strip(t):
                t = tuple(t)
                while t and t[-1] == 0:
                    t = t[:-1]
                return t
            if strip(ua) != strip(lib_alpha):
                mismatches += 1
                if mismatches <= 15:
                    print(f"MISMATCH perm={perm.trimcode} perm2={perm2.trimcode} "
                          f"lam={lam} lib_alpha={lib_alpha} user_alpha={ua}")
    print(f"\nn={n}  RC crystals checked={total}  extremal_weight != user_alpha : {mismatches}")


def strip(t):
    t = tuple(t)
    while t and t[-1] == 0:
        t = t[:-1]
    return t


def rc_anchor(hw):
    """Unique RC element of the crystal carrying the extremal weight alpha."""
    alpha = strip(hw.extremal_weight)
    for c in hw.full_crystal:
        if strip(c.length_vector) == alpha:
            return c
    return None


def wc_crystals(perm, length):
    """WC crystals of Gx(perm), grouped by WC's own highest weight; keyed by alpha."""
    by_hw = defaultdict(list)
    for wc in WCGraph.all_wc_graphs(perm, length):
        by_hw[wc.to_highest_weight()[0]].append(wc)
    out = {}
    for hw, crystal in by_hw.items():
        alpha = strip(hw.extremal_weight)
        anchor = None
        for c in hw.full_crystal:
            if strip(c.length_vector) == alpha:
                anchor = c
                break
        out.setdefault(alpha, []).append((hw, list(hw.full_crystal), anchor))
    return out


def anchor_iso(n):
    """Direct RC-anchor <-> WC-anchor bijection keyed on the extremal weight alpha."""
    length = n + 2
    total_rc = 0
    unmatched = 0
    char_mismatch = 0
    for perm in Permutation.all_permutations(n):
        if perm.inv == 0:
            continue
        # collect RC crystals keyed by alpha
        rc_by_alpha = defaultdict(list)
        for perm2, coeff, hw, full, lib_alpha, weights, lam in rc_crystal_data(perm, length):
            rc_by_alpha[strip(lib_alpha)].append((hw, full, coeff))
        wc_by_alpha = wc_crystals(perm, length)
        printed = False
        for alpha, rc_list in sorted(rc_by_alpha.items(), key=lambda t: str(t[0])):
            for hw, full, coeff in rc_list:
                total_rc += 1
                r_anchor = rc_anchor(hw)
                wc_match = wc_by_alpha.get(alpha)
                if not wc_match:
                    unmatched += 1
                    if unmatched <= 10:
                        print(f"UNMATCHED perm={perm.trimcode} alpha={alpha} (no WC crystal)")
                    continue
                # compare characters (multiset of weights, sorted) of first WC match
                wc_hw, wc_full, wc_anchor = wc_match[0]
                rc_wts = sorted(strip(c.length_vector) for c in full)
                wc_wts = sorted(strip(c.length_vector) for c in wc_full)
                iso = (rc_wts == wc_wts)
                if not iso:
                    char_mismatch += 1
                    print(f"CHAR-MISMATCH perm={perm.trimcode} alpha={alpha}  "
                          f"|RC crystal|={len(full)} |WC crystal|={len(wc_full)}  "
                          f"anchors_match={strip(r_anchor.length_vector)==strip(wc_anchor.length_vector) if (r_anchor and wc_anchor) else None}")
                if perm.inv > 0 and not printed and len(rc_by_alpha) > 0:
                    printed = True
        # per-perm summary line
    print(f"\nn={n}: RC crystals={total_rc}  unmatched_alpha={unmatched}  "
          f"char_mismatch(WC!=RC)={char_mismatch}")


if __name__ == "__main__":
    import sys
    mode = sys.argv[2] if len(sys.argv) > 2 else "verify"
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    if mode == "iso":
        anchor_iso(n)
    else:
        main(n)
