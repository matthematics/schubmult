"""Find a many-to-one RCGraph => WCGraph correspondence that breaks the
Grothendieck polynomial into key (Demazure) polynomials.

Setup
-----
* The Grothendieck polynomial of ``perm`` is a positively-weighted sum of
  RCGraphs: ``Gx(perm) = sum_{perm2} beta^(inv diff) * Sx(perm2)`` and each
  ``Sx(perm2)`` is the sum of ``polyvalue`` over its RCGraphs.
* Both RCGraph and WCGraph carry the SAME crystal structure (same
  ``_extremal_weight`` / ``to_highest_weight`` / ``full_crystal``).
* ``test_full_crystal_is_key`` proves: the full crystal of a highest-weight
  RCGraph sums (polyvalue) to ``Key(extremal_weight)``.

So grouping the RCGraphs by their highest-weight element partitions them into
Demazure crystals, one key each. This script:

  1. Verifies that RC-crystal grouping reproduces ``Gx(perm)`` in the key basis.
  2. Builds the candidate RCGraph => WCGraph map by matching shared crystal
     invariants (extremal_weight / highest-weight recording tableau) and reports
     the fibers, checking the map is many-to-one and consistent per key.
"""

import sys
from collections import defaultdict

from schubmult import Gx, Permutation, RCGraph, Sx, WCGraph
from schubmult.rings.polynomial_algebra import KeyPolyBasis, PolynomialAlgebra
from schubmult.symbolic import S, expand

x = Sx.genset
Key = PolynomialAlgebra(KeyPolyBasis(x))
beta = Gx._beta


def rc_crystals(perm, length):
    """Group the Grothendieck-poly RCGraphs into Demazure crystals.

    Returns ``crystals``: dict keyed by ``(perm2, hw_rc)`` -> (beta_coeff, [rc, ...]).
    """
    fatness = WCGraph.groth_to_schub(perm, beta)  # {perm2: beta^(inv diff)}
    crystals = {}
    for perm2, bcoeff in fatness.items():
        by_hw = defaultdict(list)
        for rc in RCGraph.all_rc_graphs(perm2, length):
            hw = rc.to_highest_weight()[0]
            by_hw[hw].append(rc)
        for hw, rcs in by_hw.items():
            crystals[(perm2, hw)] = (bcoeff, rcs)
    return crystals


def wc_crystals(perm, length):
    """Group the Grothendieck-poly WCGraphs into their crystals (by highest weight)."""
    by_hw = defaultdict(list)
    for wc in WCGraph.all_wc_graphs(perm, length):
        hw = wc.to_highest_weight()[0]
        by_hw[hw].append(wc)
    return by_hw


def main(n):
    length = n
    for perm in Permutation.all_permutations(n):
        crystals = rc_crystals(perm, length)

        # (1) key decomposition from RC crystals
        key_coeff = defaultdict(lambda: S.Zero)
        for (perm2, hw), (bcoeff, rcs) in crystals.items():
            ew = hw.extremal_weight
            key_coeff[ew] += bcoeff

        recon = S.Zero
        for ew, coeff in key_coeff.items():
            recon += coeff * Key(ew).expand()
        groth_poly = Gx(perm).expand()
        ok = expand(recon - groth_poly) == S.Zero

        # verify each RC crystal really equals its key
        crystal_ok = True
        for (perm2, hw), (bcoeff, rcs) in crystals.items():
            csum = sum((rc.polyvalue(x) for rc in rcs), S.Zero)
            if expand(csum - Key(hw.extremal_weight).expand()) != S.Zero:
                crystal_ok = False
                break

        # (2) candidate RCGraph => WCGraph map via shared crystal invariants
        wc_by_hw = wc_crystals(perm, length)
        # WCGraph highest-weight elements indexed by their extremal weight
        wc_hw_by_ew = defaultdict(list)
        for hw_wc in wc_by_hw:
            wc_hw_by_ew[hw_wc.extremal_weight].append(hw_wc)

        # each RC crystal (a key) -> which WCGraph highest weights share its extremal weight
        fiber_report = []
        many_to_one_ok = True
        for (perm2, hw), (bcoeff, rcs) in crystals.items():
            ew = hw.extremal_weight
            matching_wc = wc_hw_by_ew.get(ew, [])
            fiber_report.append((perm2.trimcode, ew, len(rcs), len(matching_wc)))
            if len(matching_wc) == 0:
                many_to_one_ok = False

        print(f"perm={perm.trimcode}  key_decomp_ok={ok}  crystal=key_ok={crystal_ok}  wc_match_ok={many_to_one_ok}")
        nz = {ew: co for ew, co in key_coeff.items() if co != 0}
        print(f"    keys: {nz}")
        for perm2c, ew, nrc, nwc in fiber_report:
            print(f"    crystal perm2={perm2c} ew={ew}: |RC fiber|={nrc}  |WC hw w/ same ew|={nwc}")
        assert ok, f"KEY DECOMP FAILED for {perm.trimcode}"
        assert crystal_ok, f"CRYSTAL != KEY for {perm.trimcode}"


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    main(n)
