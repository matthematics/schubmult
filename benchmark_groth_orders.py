"""Wall-clock benchmark of both layer orders for unpruned symbolic mixed Grothendieck products.

Example (Linux): PYTHONHASHSEED=0 OPENBLAS_NUM_THREADS=1 taskset -c 2 \
    python benchmark_groth_orders.py --repeats 7 --csv timings.csv
Use --case C:21:3142 or Q:21:3142 to select pairs; --normalize includes
the optional final coefficient normalization. No intermediate polynomial zero
tests or numerical probes are used. Imports and correctness checks are untimed.
"""

import argparse
import csv
import gc
import hashlib
import math
import multiprocessing as mp
import random
import resource
import sys
import time
from functools import cache

from schubmult import Permutation
from schubmult.abc import beta, y, z
from schubmult.combinatorics.permutation import uncode
from schubmult.mult import groth_double as gd
from schubmult.mult import groth_quantum_double as gq
from schubmult.symbolic import S
from schubmult.symbolic.poly.schub_poly import call_zvars
from schubmult.utils.schub_lib import compute_vpathdicts


@cache
def numerator(alphabet, indices, diff, power):
    """The same factorial elementary subset DP as the kernel, without probes."""
    p = len(alphabet) - diff
    state = {(0, 0): S.One}
    for i, value in enumerate(alphabet):
        nxt = {}
        for (c, r), acc in state.items():
            if i - c < diff:
                nxt[c, r] = nxt.get((c, r), S.Zero) + acc
            if c < p:
                zz = z[indices[i - c]]
                key = (c + 1, r + (value is not None))
                factor = y[value] - zz if value is not None else -(1 + beta * zz)
                nxt[key] = nxt.get(key, S.Zero) + acc * factor
        state = nxt
    return sum((beta ** (power + r) * a for (c, r), a in state.items() if c == p), S.Zero)


@cache
def edge(k, i, w, wp, v1, v2, diff, length):
    alphabet, den, left = gd._tilde_window_fate(k, w, wp)
    if diff > len(alphabet) or length < left + alphabet.count(None):
        return None  # Structural inadmissibility, not a test on a computed coefficient.
    indices = tuple(call_zvars(v1, v2, k, i)[: diff + 1])
    a = numerator(alphabet, indices, diff, length - k + diff)
    return (-a if left % 2 else a), den


def add(f, g):
    if f is None:
        return g
    a, da = f
    b, db = g
    if da == db:
        return a + b, da
    d = da.copy()
    for i, n in db.items():
        if n > d.get(i, 0):
            d[i] = n
    for i, n in d.items():
        na, nb = da.get(i, 0), db.get(i, 0)
        if n > na:
            a *= gd._atom_pow(y, beta, i, n - na)
        if n > nb:
            b *= gd._atom_pow(y, beta, i, n - nb)
    return a + b, d


def mul(f, g, weight):
    a, da = f
    b, db = g
    n = a * b * weight
    if not da:
        return n, db
    if not db:
        return n, da
    d = da.copy()
    for i, p in db.items():
        d[i] = d.get(i, 0) + p
    return n, d


def unpruned(u, v, mode, wide):
    quantum = mode == "Q"
    th = tuple(k for k in ((~v).strict_theta() if quantum else (~v).theta()) if k)
    mu = uncode(th)
    tops = {xi * mu: c for xi, c in gd.dgroth_phantom_expansion(v, z, beta, th).items()}
    graph = [{} for _ in th]
    for top in tops:
        for index, layer in enumerate(compute_vpathdicts(tuple(th), top)):
            for v1, steps in layer.items():
                for v2, d, s in steps:
                    src, dst = (v1, v2) if wide else (v2, v1)
                    graph[index].setdefault(src, set()).add((dst, v1, v2, d, int(s)))
    identity = Permutation([])
    front = {u: {identity: (S.One, {})} if wide else {t: (c, {}) for t, c in tops.items()}}
    for index in range(len(th)) if wide else reversed(range(len(th))):
        k = th[index]
        layer = graph[index]
        nxt = {}
        for w, vals in front.items():
            live = [(f, layer[node]) for node, f in vals.items() if node in layer]
            if not live:
                continue
            if quantum:
                support = ((wp, length, gq._qmon(vec, gq._vars.q_var)) for wp, (length, vec) in gq.quantum_pieri_chains(w, k).items())
            else:
                support = ((wp, wp.inv - w.inv, 1) for wp in gd._top_block_support(w, k) | {w})
            for wp, length, weight in support:
                bucket = None
                for f, steps in live:
                    for dst, v1, v2, d, s in steps:
                        c = edge(k, index + 1, w, wp, v1, v2, d, length)
                        if c is None:
                            continue
                        if bucket is None:
                            bucket = nxt.setdefault(wp, {})
                        bucket[dst] = add(bucket.get(dst), mul(f, c, s * weight))
        front = nxt
    out = {}
    for w, vals in front.items():
        if wide:
            total = None
            for t, (a, d) in vals.items():
                total = add(total, (a * tops[t], d))
            a, d = total
        else:
            a, d = vals[identity]
        out[w] = a / math.prod(gd._atom_pow(y, beta, i, p) for i, p in d.items())
    return out


def fingerprint(out):
    """Hash the support and all final expressions at two independent field points."""
    dig = hashlib.sha256()
    support = set()
    for prime, salt in [(1_000_000_007, "order-a:"), (1_000_000_009, "order-b:")]:

        @cache
        def evaluate(expr):
            if expr.is_Integer:
                return int(expr) % prime
            if expr.is_Symbol:
                return int.from_bytes(hashlib.blake2b((salt + str(expr)).encode(), digest_size=8).digest()) % prime
            if expr.is_Pow:
                a, n = expr.args
                return pow(evaluate(a), int(n), prime)
            if expr.is_Add:
                return sum(map(evaluate, expr.args)) % prime
            if expr.is_Mul:
                return math.prod(map(evaluate, expr.args)) % prime
            if expr.is_Rational:
                a, b = expr.get_num_den()
                return int(a) * pow(int(b), -1, prime) % prime
            raise ValueError(f"Unsupported expression type {type(expr)}")

        for w, c in sorted(out.items(), key=lambda item: tuple(item[0])):
            value = evaluate(c)
            if value:
                support.add(w)
                dig.update(f"{tuple(w)}:{value};".encode())
        evaluate.cache_clear()
    return len(support), dig.hexdigest()


def worker(case, order, normalize, pipe, verify):
    mode, left, right = case.split(":")
    u, v = (Permutation(list(map(int, p))) for p in (left, right))

    def product():
        res = unpruned(u, v, mode, order == "wide")
        return {w: gd.normalize_coeff(c, y, beta) for w, c in res.items()} if normalize else res

    gc.collect()
    t = time.perf_counter_ns()
    out = product()
    cold = (time.perf_counter_ns() - t) / 1e9
    del out
    loops = min(32, max(1, round(0.020 / cold)))
    t = time.perf_counter_ns()
    for _ in range(loops):
        out = product()
        del out
    warm = (time.perf_counter_ns() - t) / (1e9 * loops)
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    size, digest = fingerprint(product())
    if verify:
        func = gd.grothmult_double if mode == "C" else gq.grothmult_q_double
        if (size, digest) != fingerprint(func({u: S.One}, v, y, z, beta)):
            raise ValueError(f"Output mismatch with independent public kernel: {case} {order}")
    pipe.send((cold, warm, loops, size, digest, rss))
    pipe.close()


DEFAULT_CASES = """
C:21:3142 C:21:15243 C:21:51243 C:21:45312 C:21:54321 C:321:3142
C:31524:45312 C:35142:51243 C:24153:15243 C:21543:54321
C:21:615243 C:35142:645312 C:321:654321 C:35142:7654321
Q:21:3142 Q:321:3142 Q:21:45312 Q:21:51243 Q:31524:45312 Q:21:54321
""".split()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--case", action="append", help="C:u:v or Q:u:v (one-line inputs of rank <=9)")
    p.add_argument("--repeats", type=int, default=7)
    p.add_argument("--timeout", type=float, default=45, help="wall seconds per worker, including checks and warm runs")
    p.add_argument("--csv", required=True)
    p.add_argument("--normalize", action="store_true")
    opts = p.parse_args()
    cases = list(opts.case or DEFAULT_CASES)
    if opts.repeats < 1 or any(c.split(":")[0] not in ("C", "Q") for c in cases):
        p.error("need positive repeats and modes C or Q")
    rng = random.Random(3142)
    ctx = mp.get_context("fork")  # Parent never computes a product; each timed child starts with the same imported modules.
    blocked = set()
    expected = {}
    with open(opts.csv, "x", newline="") as dest:
        table = csv.writer(dest)
        table.writerow(["case", "repeat", "order", "normalized", "cold_s", "warm_s", "warm_loops", "support", "digest", "worker_peak_kib", "status"])
        for rep in range(opts.repeats):
            rng.shuffle(cases)
            for case in cases:
                orders = ["small", "wide"]
                if (rep + sum(map(ord, case))) % 2:
                    orders.reverse()
                for order in orders:
                    if (case, order) in blocked:
                        continue
                    reader, writer = ctx.Pipe(duplex=False)
                    proc = ctx.Process(target=worker, args=(case, order, opts.normalize, writer, rep == 0))
                    proc.start()
                    writer.close()
                    proc.join(opts.timeout)
                    if proc.is_alive():
                        proc.terminate()
                        proc.join()
                        blocked.add((case, order))
                        table.writerow([case, rep, order, opts.normalize, "", "", "", "", "", "", "timeout"])
                        print(case, order, "timeout", flush=True)
                    else:
                        if proc.exitcode != 0:
                            raise RuntimeError(f"{case} {order}: worker exit {proc.exitcode}")
                        cold, warm, loops, size, digest, rss = reader.recv()
                        if expected.setdefault(case, (size, digest)) != (size, digest):
                            raise ValueError(f"Output mismatch for {case} {order}")
                        table.writerow([case, rep, order, opts.normalize, cold, warm, loops, size, digest, rss, "ok"])
                        if not rep:
                            print(case, order, "cold/warm ms:", round(1000 * cold, 2), round(1000 * warm, 2), "outputs", size, flush=True)
                    reader.close()
                    proc.close()
                    dest.flush()
            print("completed round", rep + 1, flush=True)
    return 1 if blocked else 0


if __name__ == "__main__":
    sys.exit(main())
