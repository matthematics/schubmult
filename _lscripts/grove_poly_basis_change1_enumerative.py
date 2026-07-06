"""Enumerative check for a forest-polynomial decomposition of grove polynomials.

For each permutation w in S_n, we:
1. Build all WC graphs of w.
2. Group the forest-compatible graphs by forest weight.
3. Form the weighted forest-polynomial sum.
4. Compare that sum against the grove polynomial of w.

This script is intended as a shareable verification utility (clear output + fail-fast option).
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass

from schubmult import *
from schubmult.combinatorics.indexed_forests import *
from schubmult.combinatorics.pipe_dream import PipeDream
from schubmult.rings.combinatorial.forest_rc_ring import ForestRCGraphRing
from schubmult.rings.polynomial_algebra import *
from schubmult.symbolic.common_polys import *
from schubmult.symbolic.poly.variables import *
from schubmult.utils._mul_utils import add_perm_dict
from schubmult.utils.tuple_utils import pad_tuple

def grove_to_forest_as_rc(comp: Permutation):
    boip = WCGraph.grove_wcs(comp)
    #spits = [PipeDream.from_rc_graph(bpd).co_pipe_dream() for bpd in boip if PipeDream.from_rc_graph(bpd).co_pipe_dream().is_reduced]
    ret = {}
    for wc in boip:
        #print(wc)
        cpdb = PipeDream.from_rc_graph(wc).co_pipe_dream()
        if cpdb.is_reduced:
            w0 = Permutation.w0(cpdb.rows)
            perm = cpdb.perm * w0
            code = tuple(perm.pad_code(len(comp)))
            ret[code] = ret.get(code, 0) + 1
    return ret


@dataclass(frozen=True)
class CheckResult:
    composition: tuple[int, ...]
    wc_count: int
    forest_term_count: int
    is_match: bool
    lhs: object
    rhs: object
    diff: object


# def _forest_weight_multiplicities(wc_graphs):
#     """Count forest weights for WC graphs satisfying forest_weight == length_vector."""
#     multiplicities = {}
#     for wc_graph in wc_graphs:
#         if wc_graph.forest_weight != wc_graph.length_vector:
#             continue
#         forest_weight = wc_graph.forest_weight
#         multiplicities[forest_weight] = multiplicities.get(forest_weight, 0) + 1
#     return multiplicities


def _weighted_forest_sum(base_composition, multiplicities, beta):
    """Build the weighted forest-polynomial sum for a fixed base composition."""
    return sum(
        [
            beta ** (sum(forest_weight) - sum(base_composition))
            * multiplicity
            * Forest(*forest_weight).expand()
        
        for forest_weight, multiplicity in multiplicities.items()
        ]
    )


def verify_composition(composition, rank):
    """Verify one composition-level identity and return full diagnostic data."""
    wc_graphs = WCGraph.grove_wcs(composition)
    multiplicities = grove_to_forest_as_rc(composition)

    lhs = _weighted_forest_sum(composition, multiplicities, beta=Gx._beta)
    rhs = grove_polynomial(composition, Sx.genset, Gx._beta)
    diff = (lhs - rhs).expand()

    return CheckResult(
        composition=composition,
        wc_count=0,
        forest_term_count=len(multiplicities),
        is_match=(diff == 0),
        lhs=lhs,
        rhs=rhs,
        diff=diff,
    )


def _format_result_line(index, total, result: CheckResult):
    status = "OK" if result.is_match else "FAIL"
    return (
        f"[{index}/{total}] {status} comp={result.composition} "
        f"wc={result.wc_count} forest_terms={result.forest_term_count}"
    )


def run_checks(n, verbose=False, fail_fast=False):
    """Run the enumerative check across all permutations in S_n."""
    rank = n - 1
    permutations = Permutation.all_permutations(n)
    compositions = [perm.pad_code(rank) for perm in permutations]

    failures = []

    print(f"Checking {len(compositions)} compositions from S_{n} (rank={rank}).")
    for index, composition in enumerate(compositions, start=1):
        result = verify_composition(composition, rank)
        if verbose or not result.is_match:
            print(_format_result_line(index, len(compositions), result))

        if not result.is_match:
            failures.append(result)
            print("  LHS:", result.lhs)
            print("  RHS:", result.rhs)
            print("  DIFF:", result.diff)
            if fail_fast:
                break

    success_count = len(compositions) - len(failures)
    print(f"Summary: {success_count}/{len(compositions)} checks passed.")

    return failures


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Enumerative verification that a weighted forest-polynomial decomposition "
            "matches the Grove polynomial."
        )
    )
    parser.add_argument("n", type=int, help="Compute over permutations in S_n.")
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print one status line for every composition.",
    )
    parser.add_argument(
        "--fail-fast",
        action="store_true",
        help="Stop at the first mismatch.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    if args.n < 1:
        raise SystemExit("n must be at least 1")

    failures = run_checks(args.n, verbose=args.verbose, fail_fast=args.fail_fast)
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
