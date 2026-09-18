"""CLI: dual forest polynomial computations from arXiv:2306.10939 (Nadeau-Tewari).

Two formulas are supported, both expressed as elements of
``FreeAlgebra(WordBasis)`` keyed by composition (exponent) tuples:

  * ``--mode dual`` (default):  signed lower-ideal expansion, equal to P_F by
    Theorem 4.1::

        P_F = sum_{lower ideals L of F} (-1)^|L|
                  * sum_{L-compatible kappa} x^kappa

  * ``--mode tilde``:  the unsigned dual forest polynomial of eq. (4.1)::

        ~P_F = sum_{internal(F)-compatible kappa} x^kappa

  * ``--mode both``:  print both.
"""

import argparse
import sys

from schubmult.combinatorics.indexed_forests import (
    dual_forest_polynomial,
    minimum_n_for_dual_forest,
    tilde_forest_polynomial,
    weak_composition_to_indfor,
)
from schubmult import *

def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="dual_forest_signed_expansion",
        description=(
            "Compute the dual forest polynomial of arXiv:2306.10939 in "
            "FreeAlgebra(WordBasis). 'dual' is the signed Theorem-4.1 sum "
            "(equals P_F); 'tilde' is the unsigned eq.(4.1) polynomial ~P_F."
        ),
    )
    parser.add_argument(
        "composition",
        nargs="+",
        type=int,
        help="Weak composition entries (space-separated), e.g. 2 0 1",
    )
    parser.add_argument(
        "-n",
        "--num-variables",
        type=int,
        dest="n",
        help="Restrict labelings to [1..n]; auto-picked if omitted",
    )
    parser.add_argument(
        "--mode",
        choices=["dual", "tilde", "both"],
        default="dual",
        help="dual=Thm 4.1 signed sum (default), tilde=eq.(4.1), both=print both",
    )
    parser.add_argument(
        "--raw-dict",
        action="store_true",
        help="Also print the underlying coefficient dict",
    )

    args = parser.parse_args(argv)
    if args.n is not None and args.n <= 0:
        raise ValueError("--num-variables must be positive")

    code = tuple(args.composition)
    forest = weak_composition_to_indfor(code)

    n_min = minimum_n_for_dual_forest(forest)
    n = args.n if args.n is not None else n_min
    if args.mode in {"tilde", "both"} and n < n_min:
        print(f"Note: n={n} too small for ~P_F; using n={n_min} for tilde mode.")
    n_tilde = max(n, n_min) if args.mode in {"tilde", "both"} else n

    print(f"Composition: {code}    n = {n}")
    print(ForestDual(*code).change_basis(WordBasis))
    if args.mode in {"dual", "both"}:
        elem = dual_forest_polynomial(forest, n)
        print("Dual forest polynomial P_F (Thm 4.1, signed sum) in FreeAlgebra(WordBasis):")
        print(elem)
        
        if args.raw_dict:
            print("dual dict:", dict(sorted(elem.items())))

    if args.mode in {"tilde", "both"}:
        elem = tilde_forest_polynomial(forest, n_tilde)
        print(f"~P_F (eq. 4.1, unsigned, n = {n_tilde}) in FreeAlgebra(WordBasis):")
        print(elem)
        print(elem.change_basis(ForestBasis))
        if args.raw_dict:
            print("tilde dict:", dict(sorted(elem.items())))


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
