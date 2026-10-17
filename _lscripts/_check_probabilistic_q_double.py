"""Check schubmult_q_double_fast(probabilistic=True) against the exact kernel: S_4 x S_4 exhaustively by
expansion (C++ and Python kernels, same and mixed alphabets), then timing and agreement on larger products
via the probabilistic output versus the exact output filtered by expansion."""

import itertools
import resource
import time

resource.setrlimit(resource.RLIMIT_AS, (6 << 30, 6 << 30))

from schubmult.combinatorics.permutation import Permutation  # noqa: E402
from schubmult.mult.quantum_double import _schubmult_q_double_fast_python, schubmult_q_double_fast  # noqa: E402
from schubmult.symbolic import expand, sympify  # noqa: E402
from schubmult.symbolic.common_polys import _vars  # noqa: E402
from schubmult.symbolic.poly.variables import GeneratingSet  # noqa: E402

y, z, q = GeneratingSet("y"), GeneratingSet("z"), _vars.q_var


def by_expansion(d):
    return {w: expand(c) for w, c in d.items() if expand(c) != 0}


S4 = [Permutation(p) for p in itertools.permutations([1, 2, 3, 4])]
bad = 0
t = time.perf_counter()
for a, b in itertools.product(S4, S4):
    for v3 in (y, z):
        e = by_expansion(schubmult_q_double_fast({a: sympify(1)}, b, y, v3, q))
        p1 = by_expansion(schubmult_q_double_fast({a: sympify(1)}, b, y, v3, q, probabilistic=True))
        p2 = by_expansion(_schubmult_q_double_fast_python({a: sympify(1)}, b, y, v3, q, probabilistic=True))
        bad += (e != p1) + (e != p2)
print(f"S4 x S4 mismatches: {bad} ({time.perf_counter() - t:.1f}s)", flush=True)

for u, v in (([3, 1, 5, 2, 4], [2, 5, 1, 4, 3]), ([4, 1, 6, 5, 2, 3], [3, 6, 1, 5, 2, 4]), ([4, 1, 6, 5, 2, 3], [6, 1, 5, 4, 2, 3])):
    for v3, label in ((y, "same"), (z, "mixed")):
        t = time.perf_counter()
        exact = schubmult_q_double_fast({Permutation(u): sympify(1)}, Permutation(v), y, v3, q)
        t_exact = time.perf_counter() - t
        t = time.perf_counter()
        prob = schubmult_q_double_fast({Permutation(u): sympify(1)}, Permutation(v), y, v3, q, probabilistic=True)
        t_prob = time.perf_counter() - t
        t = time.perf_counter()
        ok = by_expansion(exact) == by_expansion(prob)
        print(f"{u} x {v} ({label}): exact {t_exact:.2f}s ({len(exact)} terms), probabilistic {t_prob:.2f}s ({len(prob)} terms), equal after expansion: {ok} ({time.perf_counter() - t:.1f}s)", flush=True)
