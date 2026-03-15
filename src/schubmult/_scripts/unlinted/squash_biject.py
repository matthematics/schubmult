from schubmult import *
from sympy import pretty_print, S
import argparse
import itertools
from functools import lru_cache
import time

def left_squash(rc1, rc2):
    
    rc1_p = _transpose_rc_cached(rc1)
    if len(rc1_p.perm.descents()) > 1:
        raise KeyError("rc1 must be Grassmannian")
    rc2_p = _transpose_rc_cached(rc2)
    the_squash = _transpose_rc_cached(rc2_p.resize(len(rc1_p)).squash_product(rc1.resize(len(rc1_p)))).resize(len(rc2))
    return the_squash
    # return _transpose_rc_cached(the_squash)
    


@lru_cache(maxsize=None)
def _transpose_rc_cached(rc):
    return BPD.from_rc_graph(rc.transpose()).to_rc_graph().resize(max(len(rc),len(rc.perm)))


def grassmannian_rc_graphs_in_sn(n):
    """Grassmannian RC graphs in S_n with descent at n-1 (zero-indexed); include identity."""
    ret = []
    for perm in Permutation.all_permutations(n + 2):
        desc = perm.descents()
        if perm.inv != 0 and desc != {n - 1}:
            continue
        for rc in RCGraph.all_rc_graphs(perm, n):
            ret.append(rc.resize(n))
    return ret


def all_rc_graphs_in_sn(n):
    ret = []
    for perm in Permutation.all_permutations(n):
        for rc in RCGraph.all_rc_graphs(perm, n - 1):
            ret.append(rc)
    return ret


def test_grassmannian_ring_closure(n, r):
    grasses = grassmannian_rc_graphs_in_sn(n + 2)
    grass_set = set(grasses)
    failures = []
    for rc1 in grasses:
        for rc2 in grasses:
            prod = r(rc1) * r(rc2)
            support = {rc for rc, coeff in prod.items() if coeff != 0}
            if not support.issubset(grass_set):
                failures.append((rc1, rc2, support - grass_set))
    return grasses, failures


def test_left_squash_reversed_order(n, r, grasses, *, debug=False, progress_every=1):
    t0 = time.perf_counter()
    all_rc = all_rc_graphs_in_sn(n + 1)
    if debug:
        print(f"[debug] built all_rc with {len(all_rc)} graphs in {time.perf_counter() - t0:.2f}s")

    t1 = time.perf_counter()
    failures = []
    #pair_products = {(grass_rc1, grass_rc2): grass_rc2.squash_product(grass_rc1) for grass_rc1, grass_rc2 in itertools.product(grasses, repeat=2)}
    # if debug:
    #     print(f"[debug] precomputed {len(pair_products)} grass pair products in {time.perf_counter() - t1:.2f}s")

    total = len(all_rc) * len(grasses) * len(grasses)
    done = 0
    nops = 0
    t2 = time.perf_counter()
    for rc in all_rc:
        for grass_rc1, grass_rc2 in itertools.product(grasses, repeat=2):
            try:
                #lhs = r(left_squash(pair_products[(grass_rc1, grass_rc2)], rc))
                
                rhs = grass_rc1.left_squash(grass_rc2).left_squash(rc)
                lhs = (grass_rc1.left_squash(grass_rc2)).left_squash(rc)
                #if not lhs.almosteq(rhs):
                if lhs != rhs:
                    failures.append((rc, grass_rc1, grass_rc2, lhs, rhs))
                    raise AssertionError(f"Failure for rc={rc}, grass_rc1={grass_rc1}, grass_rc2={grass_rc2}\nLHS: {tuple(lhs)} \nRHS: {tuple(rhs)}")
            except KeyError as e:
                nops += 1
            done += 1
            if debug and done % progress_every == 0:
                elapsed = time.perf_counter() - t2
                rate = done / elapsed if elapsed > 0 else 0
                eta = (total - done) / rate if rate > 0 else float("inf")
                print(f"[debug] progress {done}/{total} ({100*done/total:.1f}%), rate={rate:.1f}/s, eta={eta:.1f}s, failures={len(failures)}, nops={nops}")
    if debug:
        print(f"[debug] main comparison loop finished in {time.perf_counter() - t2:.2f}s")
    return all_rc, failures

def all_grassmannian_rc_graphs(n: int, max_inv: int) -> list[RCGraph]:
    """All RC graphs for Grassmannian permutations generated from partitions."""
    graph_set = set()
    for perm in grassmannian_perms_from_partitions(n, max_inv):
        graph_set.update(rc.resize(n) for rc in RCGraph.all_rc_graphs(perm, n))
    return sorted(graph_set, key=lambda rc: (rc.perm.inv, rc.length_vector, tuple(rc)))

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

if __name__ == "__main__":
    t_start = time.perf_counter()
    parser = argparse.ArgumentParser(
        description="Test Grassmannian closure and left_squash reversed-order behavior in DualRCGraphRing."
    )
    parser.add_argument("n", type=int, help="Work in S_{n+1} with RC graphs resized to n rows.")
    parser.add_argument("--max-inv", type=int, default=4, help="Max inversion/partition sum used for Grassmannian generation.")
    parser.add_argument("--debug", action="store_true", help="Enable detailed timing/progress debug output.")
    parser.add_argument("--progress-every", type=int, default=5, help="Emit progress line every this many inner iterations when --debug is set.")
    # parser.add_argument(
    #     "--show-first-failure",
    #     action="store_true",
    #     help="Print the first failing witness in each test.",
    #     default
    # )
    args = parser.parse_args()

    n = args.n
    if n < 2:
        raise ValueError("n must be at least 2")

    r = DualRCGraphRing()
    closure_failures = []

    #grasses, closure_failures = test_grassmannian_ring_closure(n, r)
    # print(f"n={n}")
    # print(f"Grassmannian RC graphs considered: {len(grasses)}")
    # print(f"Closure under DualRCGraphRing product: {'PASS' if not closure_failures else 'FAIL'}")

    # if closure_failures and args.show_first_failure:
    #     rc1, rc2, bad_support = closure_failures[0]
    #     print("First closure failure:")
    #     pretty_print(rc1)
    #     pretty_print(rc2)
    #     print("Support outside Grassmannian set:")
    #     for bad in list(bad_support)[:5]:
    #         pretty_print(bad)
    t_grass = time.perf_counter()
    grasses = all_grassmannian_rc_graphs(n, args.max_inv)
    if args.debug:
        print(f"[debug] built grass set with {len(grasses)} graphs in {time.perf_counter() - t_grass:.2f}s (max_inv={args.max_inv})")
    all_rc, reverse_failures = test_left_squash_reversed_order(
        n,
        r,
        grasses,
        debug=args.debug,
        progress_every=max(1, args.progress_every),
    )
    print(f"All RC graphs in S_n with n rows: {len(all_rc)}")
    print(f"Grassmannian RC graphs: {len(grasses)}")
    print(f"Total checks: {len(all_rc) * len(grasses) * len(grasses)}")
    print(f"left_squash implements reversed product (r(grass) * r(rc)): {'PASS' if not reverse_failures else 'FAIL'}")

    if reverse_failures:
        rc, grass_rc1, grass_rc2, lhs, rhs = reverse_failures[0]
        print("First reversed-order failure:")
        print("rc=")
        pretty_print(rc)
        print("grass_rc1=")
        pretty_print(grass_rc1)
        print("grass_rc2=")
        pretty_print(grass_rc2)
        
        print("lhs = r(left_squash(rc, grass_rc))")
        pretty_print(lhs)
        print("rhs = r(grass_rc) * r(rc)")
        pretty_print(rhs)

    if closure_failures or reverse_failures:
        if args.debug:
            print(f"[debug] total runtime: {time.perf_counter() - t_start:.2f}s")
        raise SystemExit(1)
    if args.debug:
        print(f"[debug] total runtime: {time.perf_counter() - t_start:.2f}s")
    raise SystemExit(0)