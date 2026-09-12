"""Cross-check cpp/schubmult_core against schubmult.mult.single.schubmult_py.

usage: python compare.py [n] [num_random] [seed]
Runs all pairs in S_n for n <= 5 (else random pairs), reports mismatches and timing.
"""

import itertools
import random
import subprocess
import sys
import time
from pathlib import Path

from schubmult import Permutation, schubmult_py

BIN = Path(__file__).with_name("build") / "schubmult_core"


def run_cpp(u, v):
    args = [str(BIN), *map(str, u), "-", *map(str, v)]
    out = subprocess.run(args, capture_output=True, text=True, check=True).stdout
    res = {}
    for line in out.splitlines():
        val, perm = line.split("  ", 1)
        res[Permutation(eval(perm))] = int(val)
    return res


def run_py(u, v):
    perms = sorted([Permutation(u), Permutation(v)], reverse=True, key=lambda x: sum((~x).theta()) - x.inv)
    d = schubmult_py({perms[0]: 1}, perms[1])
    return {k: val for k, val in d.items() if val != 0}


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    num_random = int(sys.argv[2]) if len(sys.argv) > 2 else 50
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    rng = random.Random(seed)
    if n <= 5:
        allp = [list(p) for p in itertools.permutations(range(1, n + 1))]
        pairs = [(u, v) for u in allp for v in allp]
    else:
        pairs = []
        for _ in range(num_random):
            u = list(range(1, n + 1))
            v = list(range(1, n + 1))
            rng.shuffle(u)
            rng.shuffle(v)
            pairs.append((u, v))
    bad = 0
    t_cpp = t_py = 0.0
    for u, v in pairs:
        t0 = time.perf_counter()
        a = run_cpp(u, v)
        t1 = time.perf_counter()
        b = run_py(u, v)
        t2 = time.perf_counter()
        t_cpp += t1 - t0
        t_py += t2 - t1
        if a != b:
            bad += 1
            print(f"MISMATCH u={u} v={v}")
            for k in set(a) | set(b):
                if a.get(k) != b.get(k):
                    print(f"   {k}: cpp={a.get(k)} py={b.get(k)}")
    print(f"{len(pairs)} pairs, {bad} mismatches; cpp {t_cpp:.2f}s (incl. process spawn), py {t_py:.2f}s")


if __name__ == "__main__":
    main()
