"""schubmult side of the exhaustive S_n timing benchmark (same record format as harness.c).

Usage: python bench_schubmult.py n outfile [workers]
Times the compiled kernel ``schubmult_cpp.schubmult_py({u: 1}, v)`` for every pair u <= v (lexicographic
index) in S_n, in ``workers`` processes sharded by u index, and writes 24-byte records
``<IIIIQ`` = (ui, vi, time_ns, nterms, hash) with hash = sum fnv1a64(trimmed w) * coeff mod 2**64.
"""
import itertools
import os
import struct
import sys
import time
from multiprocessing import Pool

MASK = (1 << 64) - 1


def fnv_perm(w):
    h = 1469598103934665603
    for a in w:
        h ^= a
        h = (h * 1099511628211) & MASK
    return h


def shard(args):
    n, start, end, path = args
    from schubmult import Permutation
    from schubmult.mult._accel import _cpp

    perms = [Permutation(list(p)) for p in itertools.permutations(range(1, n + 1))]
    hashes = [fnv_perm(tuple(p)) for p in perms]
    total = len(perms)
    out = open(path, "wb")
    pack = struct.Struct("<IIIIQ").pack
    perf = time.perf_counter_ns
    for ui in range(start, min(end, total)):
        u = perms[ui]
        for vi in range(ui, total):
            v = perms[vi]
            t0 = perf()
            res = _cpp.schubmult_py({u: 1}, v)
            t1 = perf()
            h = 0
            for w, c in res.items():
                h = (h + fnv_perm(tuple(w)) * c) & MASK
            out.write(pack(ui, vi, min(t1 - t0, 0xFFFFFFFF), len(res), h))
    out.close()
    return path


def main():
    n = int(sys.argv[1])
    outfile = sys.argv[2]
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else os.cpu_count() - 2
    total = 1
    for i in range(2, n + 1):
        total *= i
    # work for u index i is proportional to total - i: balance by interleaving small chunks
    chunk = max(1, total // (workers * 16))
    tasks = [(n, s, min(s + chunk, total), f"{outfile}.part{s:08d}") for s in range(0, total, chunk)]
    t0 = time.perf_counter()
    with Pool(workers) as pool:
        parts = pool.map(shard, tasks, chunksize=1)
    with open(outfile, "wb") as out:
        for p in sorted(parts):
            with open(p, "rb") as f:
                out.write(f.read())
            os.remove(p)
    print(f"S_{n}: {total * (total + 1) // 2} products in {time.perf_counter() - t0:.1f}s wall with {workers} workers -> {outfile}")


if __name__ == "__main__":
    main()
