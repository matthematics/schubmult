"""Check probabilistic=True against the exact kernel.

S_4 x S_4 exhaustively by expansion (C++ and Python kernels, same and mixed alphabets); the large
product by comparing the probabilistic output with the exact Sage ring in libsingular normal form
(the exact kernel output's junk coefficients must never be expanded with symengine: that exhausts
memory).
"""

import itertools
import resource
import time

resource.setrlimit(resource.RLIMIT_AS, (6 << 30, 6 << 30))  # hard cap: 6 GB

import sage.all  # noqa: E402
from sage.all import ZZ  # noqa: E402

from schubmult.combinatorics.permutation import Permutation  # noqa: E402
from schubmult.mult.double import _schubmult_double_python, schubmult_double  # noqa: E402
from schubmult.sage import DoubleSchubertPolynomialRing  # noqa: E402
from schubmult.symbolic import expand, sympify  # noqa: E402
from schubmult.symbolic.poly.variables import GeneratingSet  # noqa: E402

y, z = GeneratingSet("y"), GeneratingSet("z")


def by_expansion(d):
    return {w: expand(c) for w, c in d.items() if expand(c) != 0}


S4 = [Permutation(p) for p in itertools.permutations([1, 2, 3, 4])]
bad = 0
t = time.perf_counter()
for a, b in itertools.product(S4, S4):
    e = by_expansion(schubmult_double({a: sympify(1)}, b, y, y))
    e3 = by_expansion(schubmult_double({a: sympify(1)}, b, y, z))
    p1 = by_expansion(schubmult_double({a: sympify(1)}, b, y, y, probabilistic=True))
    p2 = by_expansion(_schubmult_double_python({a: sympify(1)}, b, y, y, probabilistic=True))
    p3 = by_expansion(_schubmult_double_python({a: sympify(1)}, b, y, z, probabilistic=True))
    p4 = by_expansion(schubmult_double({a: sympify(1)}, b, y, z, probabilistic=True))
    bad += (e != p1) + (e != p2) + (e3 != p3) + (e3 != p4)
print(f"S4 x S4 mismatches: {bad} ({time.perf_counter() - t:.1f}s)", flush=True)

u, v = [4, 1, 6, 5, 2, 3], [8, 1, 7, 6, 2, 3, 5, 4]
X = DoubleSchubertPolynomialRing(ZZ)
t = time.perf_counter()
exact = X(u) * X(v)
print(f"exact Sage ring: {time.perf_counter() - t:.2f}s, {len(exact)} terms", flush=True)
for probabilistic in (False, True):
    P = DoubleSchubertPolynomialRing(ZZ, raw_coefficients=True, probabilistic=probabilistic)
    t = time.perf_counter()
    f = P(u) * P(v)
    t_mul = time.perf_counter() - t
    t = time.perf_counter()
    ok = X(f) == exact
    print(f"raw ring probabilistic={probabilistic}: product {t_mul:.2f}s, {len(f)} terms, equals exact: {ok} (check {time.perf_counter() - t:.2f}s)", flush=True)
XP = DoubleSchubertPolynomialRing(ZZ, probabilistic=True)
t = time.perf_counter()
g = XP(u) * XP(v)
print(f"exact-coefficient ring with probabilistic kernel: {time.perf_counter() - t:.2f}s, {len(g)} terms, equals exact: {X(g) == exact}")
