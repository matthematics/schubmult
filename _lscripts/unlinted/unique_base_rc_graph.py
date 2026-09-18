from __future__ import annotations

import argparse
import json
from collections import defaultdict
from itertools import combinations

from schubmult import Permutation, RCGraph
from schubmult.combinatorics.crystal_graph import CrystalGraphTensor
from schubmult.combinatorics.permutation import uncode


def partitions_with_sum_at_most(max_sum: int, max_parts: int) -> list[tuple[int, ...]]:
    """All partitions (weakly decreasing tuples) with at most max_parts parts and total <= max_sum."""

    def rec(remaining_sum: int, max_next: int, parts_left: int):
        yield ()
        if parts_left == 0:
            return
        for first in range(1, min(remaining_sum, max_next) + 1):
            for tail in rec(remaining_sum - first, first, parts_left - 1):
                yield (first, *tail)

    seen = set()
    for part in rec(max_sum, max_sum, max_parts):
        seen.add(part)
    return sorted(seen, key=lambda p: (sum(p), len(p), p), reverse=False)


def grassmannian_perms_from_partitions(n: int, max_inv: int) -> list[Permutation]:
    """Build Grassmannian permutations from partitions: pad to length n, reverse, then uncode."""
    perms = []
    seen = set()
    for part in partitions_with_sum_at_most(max_inv, n):
        padded = (*part, *([0] * (n - len(part))))
        weakly_increasing_code = tuple(reversed(padded))
        perm = uncode(list(weakly_increasing_code))
        if perm not in seen:
            assert perm.inv == 0 or perm.descents() == set([n - 1])
            seen.add(perm)
            perms.append(perm)
    return perms


def all_grassmannian_rc_graphs(n: int, max_inv: int) -> list[RCGraph]:
    """All RC graphs for Grassmannian permutations generated from partitions."""
    graphs = []
    for perm in grassmannian_perms_from_partitions(n, max_inv):
        graphs.extend(sorted(RCGraph.all_rc_graphs(perm, n), key=lambda rc: (rc.perm.inv, rc.length_vector, tuple(rc))))
    return graphs


def all_base_rc_graphs_in_sn(n: int) -> list[RCGraph]:
    """All RC graphs for permutations in S_n, realized with n rows for right multiplication."""
    graphs = []
    for perm in Permutation.all_permutations(n):
        graphs.extend(rc.resize(n) for rc in RCGraph.all_rc_graphs(perm, n - 1))
    return sorted(graphs, key=lambda rc: (rc.perm.inv, rc.length_vector, tuple(rc)))


def tensor_signature(tensor_hw: CrystalGraphTensor) -> tuple[str, str]:
    """Canonical hashable representation of highest-weight tensor factors."""
    return (str(tensor_hw.factors[0]), str(tensor_hw.factors[1]))


def check_hw_tensor_set_disjointness(n: int, max_inv: int):
    """Group highest-weight tensors by highest-weight base RC, then test pairwise disjointness of the groups."""
    grass_graphs = all_grassmannian_rc_graphs(n, max_inv)
    base_graphs = all_base_rc_graphs_in_sn(n)

    tensors_by_base_hw: dict[RCGraph, set[tuple[str, str]]] = defaultdict(set)
    witness_by_base_hw: dict[RCGraph, dict[tuple[str, str], dict[str, str]]] = defaultdict(dict)

    for base_rc in base_graphs:
        base_hw = base_rc.to_highest_weight()[0].resize(n)
        for grass_rc in grass_graphs:
            tensor = CrystalGraphTensor(base_rc, grass_rc)
            tensor_hw, _ = tensor.to_highest_weight()
            sig = tensor_signature(tensor_hw)
            tensors_by_base_hw[base_hw].add(sig)
            if sig not in witness_by_base_hw[base_hw]:
                witness_by_base_hw[base_hw][sig] = {
                    "base_rc": str(base_rc),
                    "grass_rc": str(grass_rc),
                    "tensor_hw_left": sig[0],
                    "tensor_hw_right": sig[1],
                }

    base_hw_keys = sorted(tensors_by_base_hw.keys(), key=lambda rc: (rc.perm.inv, rc.length_vector, tuple(rc)))
    intersections = []
    for rc1, rc2 in combinations(base_hw_keys, 2):
        inter = tensors_by_base_hw[rc1].intersection(tensors_by_base_hw[rc2])
        if inter:
            sample_sig = sorted(inter)[0]
            intersections.append(
                {
                    "base_hw_1": str(rc1),
                    "base_hw_2": str(rc2),
                    "intersection_size": len(inter),
                    "sample_tensor_hw": {
                        "left": sample_sig[0],
                        "right": sample_sig[1],
                    },
                    "witness_1": witness_by_base_hw[rc1][sample_sig],
                    "witness_2": witness_by_base_hw[rc2][sample_sig],
                }
            )

    return base_hw_keys, grass_graphs, tensors_by_base_hw, intersections


def make_json_summary(base_hw_keys, grass_graphs, tensors_by_base_hw, intersections):
    return {
        "num_base_hw_keys": len(base_hw_keys),
        "num_grassmannian_rc_graphs": len(grass_graphs),
        "tensor_set_sizes": [{"base_hw": str(base_hw), "size": len(tensors_by_base_hw[base_hw])} for base_hw in base_hw_keys],
        "is_pairwise_disjoint": len(intersections) == 0,
        "intersections": [
            {
                "base_hw_1": item["base_hw_1"],
                "base_hw_2": item["base_hw_2"],
                "intersection_size": item["intersection_size"],
                "sample_tensor_hw": item["sample_tensor_hw"],
                "witness_1": item["witness_1"],
                "witness_2": item["witness_2"],
            }
            for item in intersections
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Compute highest-weight tensors for all base_rc ⊗ grass_rc, bucket by highest-weight(base_rc), "
            "and test pairwise disjointness of these tensor sets."
        )
    )
    parser.add_argument("n", type=int, help="Size n for base permutations in S_n.")
    parser.add_argument(
        "--max-inv",
        type=int,
        required=True,
        help="Maximum partition sum (equivalently inversion bound for generated Grassmannians).",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Optional JSON file to save summary data.",
    )
    args = parser.parse_args()

    if args.n < 1:
        parser.error("n must be at least 1")
    if args.max_inv < 0:
        parser.error("--max-inv must be nonnegative")

    base_hw_keys, grass_graphs, tensors_by_base_hw, intersections = check_hw_tensor_set_disjointness(args.n, args.max_inv)

    print(f"n={args.n}")
    print(f"Grassmannian RC graphs generated: {len(grass_graphs)}")
    print(f"Distinct highest-weight base_rc keys: {len(base_hw_keys)}")
    print(f"Total (base_rc, grass_rc) pairs checked: {len(all_base_rc_graphs_in_sn(args.n)) * len(grass_graphs)}")
    print("Distinct highest-weight tensors per base_hw:")
    for base_hw in base_hw_keys:
        print(f"  size {len(tensors_by_base_hw[base_hw])} for base_hw={base_hw}")

    if not intersections:
        print("True: tensor sets by highest-weight(base_rc) are pairwise disjoint.")
    else:
        print(f"False: found {len(intersections)} nonempty pairwise intersections.")
        for item in intersections:
            print(f"  intersection size {item['intersection_size']} between")
            print(f"    base_hw_1={item['base_hw_1']}")
            print(f"    base_hw_2={item['base_hw_2']}")
            print(f"    sample tensor_hw=({item['sample_tensor_hw']['left']}, {item['sample_tensor_hw']['right']})")

    if args.output is not None:
        payload = make_json_summary(base_hw_keys, grass_graphs, tensors_by_base_hw, intersections)
        with open(args.output, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)
        print(f"Wrote summary to {args.output}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
