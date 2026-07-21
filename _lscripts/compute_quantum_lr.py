from __future__ import annotations

import argparse
import itertools
import json
import platform
import random
import sys
import time
from pathlib import Path

from schubmult import *
from schubmult.rings.schubert import QDSx


def perm_to_list(perm: Permutation) -> list[int]:
    return list(perm)


def elem_to_json_terms(elem) -> list[dict[str, object]]:
    terms: list[dict[str, object]] = []
    for perm, coeff in sorted(elem.items(), key=lambda kv: tuple(perm_to_list(kv[0]))):
        terms.append({"perm": perm_to_list(perm), "coeff": str(coeff)})
    return terms


def elems_equal(elem1, elem2) -> bool:
    if hasattr(elem1, "almosteq"):
        return bool(elem1.almosteq(elem2))
    return elem1 == elem2


def run_unit_checks(perms: list[Permutation]) -> dict[str, object]:
    start = time.perf_counter()
    failures: list[dict[str, object]] = []
    schubert_identity_mismatches: list[list[int]] = []
    ring_one = QDSx(perms[0]).ring.one
    schubert_identity = QDSx(Permutation([]))
    for perm in perms:
        elem = QDSx(perm)
        left_ok = elems_equal(ring_one * elem, elem)
        right_ok = elems_equal(elem * ring_one, elem)
        if not (left_ok and right_ok):
            failures.append({
                "perm": perm_to_list(perm),
                "left_ok": left_ok,
                "right_ok": right_ok,
            })
        # This is a useful algebraic diagnostic, but not a unit-axiom failure.
        if not elems_equal(schubert_identity * elem, elem) or not elems_equal(elem * schubert_identity, elem):
            schubert_identity_mismatches.append(perm_to_list(perm))
    return {
        "checks": len(perms),
        "failed": len(failures),
        "runtime_seconds": time.perf_counter() - start,
        "failures": failures,
        "schubert_identity_mismatch_count": len(schubert_identity_mismatches),
        "schubert_identity_mismatches": schubert_identity_mismatches,
    }


def run_associativity_checks(perms: list[Permutation], samples: int, seed: int) -> dict[str, object]:
    start = time.perf_counter()
    rng = random.Random(seed)
    failures: list[dict[str, object]] = []

    for _ in range(samples):
        a, b, c = (rng.choice(perms), rng.choice(perms), rng.choice(perms))
        ea, eb, ec = QDSx(a), QDSx(b), QDSx(c)
        left = (ea * eb) * ec
        right = ea * (eb * ec)
        if not elems_equal(left, right):
            failures.append({
                "a": perm_to_list(a),
                "b": perm_to_list(b),
                "c": perm_to_list(c),
            })

    return {
        "samples": samples,
        "failed": len(failures),
        "seed": seed,
        "runtime_seconds": time.perf_counter() - start,
        "failures": failures,
    }


def compute_all_products(n: int, out_path: Path, progress_every: int) -> dict[str, object]:
    perms = Permutation.all_permutations(n)
    total_pairs = len(perms) ** 2
    pair_count = 0
    nonzero_term_count = 0
    max_terms_in_product = 0
    start = time.perf_counter()

    with out_path.open("w", encoding="utf-8") as f:
        for perm1, perm2 in itertools.product(perms, repeat=2):
            prod = QDSx(perm1) * QDSx(perm2)
            terms = elem_to_json_terms(prod)
            rec = {
                "u": perm_to_list(perm1),
                "v": perm_to_list(perm2),
                "terms": terms,
            }
            f.write(json.dumps(rec, ensure_ascii=True) + "\n")

            pair_count += 1
            num_terms = len(terms)
            nonzero_term_count += num_terms
            if num_terms > max_terms_in_product:
                max_terms_in_product = num_terms

            if progress_every > 0 and pair_count % progress_every == 0:
                elapsed = time.perf_counter() - start
                print(f"Progress: {pair_count}/{total_pairs} products in {elapsed:.1f}s", flush=True)

    runtime = time.perf_counter() - start
    return {
        "n": n,
        "permutations": len(perms),
        "total_products": total_pairs,
        "total_nonzero_terms": nonzero_term_count,
        "max_terms_in_product": max_terms_in_product,
        "runtime_seconds": runtime,
        "avg_ms_per_product": (runtime * 1000.0 / total_pairs) if total_pairs else 0.0,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compute all quantum Schubert products for FL_n and export machine-readable tables.",
    )
    parser.add_argument("n", type=int, help="Flag size n for FL_n")
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("logs/quantum_lr"),
        help="Output directory for product tables and metadata",
    )
    parser.add_argument(
        "--basename",
        default=None,
        help="Optional output base name (default: fl{n}_quantum_lr)",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=200,
        help="Print progress every k products (0 disables progress prints)",
    )
    parser.add_argument(
        "--assoc-samples",
        type=int,
        default=200,
        help="Number of random associativity checks",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Random seed for associativity sampling",
    )
    parser.add_argument(
        "--skip-checks",
        action="store_true",
        help="Skip unit and associativity checks",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.n < 1:
        parser.error("n must be >= 1")

    basename = args.basename or f"fl{args.n}_quantum_lr"
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    products_path = out_dir / f"{basename}.products.jsonl"
    metadata_path = out_dir / f"{basename}.metadata.json"

    print(f"Starting full product computation for FL_{args.n}")
    products_summary = compute_all_products(args.n, products_path, args.progress_every)

    checks: dict[str, object] = {}
    if not args.skip_checks:
        perms = Permutation.all_permutations(args.n)
        print("Running unit checks...")
        checks["unit"] = run_unit_checks(perms)
        print("Running associativity checks...")
        checks["associativity"] = run_associativity_checks(perms, args.assoc_samples, args.seed)

    metadata = {
        "version": 1,
        "command": " ".join(sys.argv),
        "python": sys.version,
        "platform": platform.platform(),
        "products_file": str(products_path),
        "products_summary": products_summary,
        "checks": checks,
    }
    metadata_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")

    print(f"Wrote products: {products_path}")
    print(f"Wrote metadata: {metadata_path}")
    if checks:
        unit_failed = checks["unit"]["failed"]
        assoc_failed = checks["associativity"]["failed"]
        print(f"Validation summary: unit_failed={unit_failed}, associativity_failed={assoc_failed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())