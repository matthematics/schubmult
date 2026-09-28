"""Benchmark schubmult's Grothendieck products against Buch's Equivariant Schubert Calculator (Maple).

Conventions (verified on Fl(3)/Fl(4)):
  X[w]  <->  G_w at beta = -1, same permutation
  T[i]  <->  (1 - y_i) / (1 - y_{i+1})
  Fl(n) drops terms with w not in S_n.

Usage:  conda activate schubmult_312; python benchmark_equivcalc.py [n] [pairs] [--equivariant]
"""
import random
import re
import subprocess
import sys
import time

import sympy

from schubmult import DGx, Gx, Permutation
from schubmult.abc import beta, y

MAPLE = "/mnt/c/Program Files/Maple 2026/bin.X86_64_WINDOWS/cmaple.exe"
EQUIVCALC = "C:/Users/matth/equivcalc/EquivCalc-1.0.1/equivcalc"
SCRIPT = "/mnt/c/Users/matth/equivcalc/bench_tmp.mpl"


def redword(perm):
    """A reduced word of perm (one-line, 1-based) as product s_{a1} s_{a2} ... left to right."""
    w = list(perm)
    word = []
    # bubble sort from the right: w = s_{a1}...s_{ak} means w(i) = s_{a1}(...(s_{ak}(i)))
    while True:
        for i in range(len(w) - 1):
            if w[i] > w[i + 1]:
                w[i], w[i + 1] = w[i + 1], w[i]
                word.append(i + 1)
                break
        else:
            break
    return word[::-1]


def maple_script(n, pairs, equivariant):
    fn = "ktmult" if equivariant else "kmult"
    lines = [f'read "{EQUIVCALC}":', "with(equivcalc):", f"Fl({n}):"]
    for k, (u, v) in enumerate(pairs):
        lines += [
            f"u := redexp_weyl({redword(u)}):",
            f"v := redexp_weyl({redword(v)}):",
            f"st := time(): r := {fn}(u, v): et := time() - st:",
            f'printf("PAIR {k} %a\\n", et):',
            'printf("RESULT %a\\n", weyl_sp(r)):',
        ]
    lines.append("quit;")
    return "\n".join(lines) + "\n"


def run_maple(script_text):
    with open(SCRIPT, "w") as f:
        f.write(script_text)
    with open(SCRIPT) as f:
        out = subprocess.run([MAPLE, "-q"], stdin=f, capture_output=True, text=True, timeout=36000).stdout
    times, results = {}, {}
    for line in out.splitlines():
        m = re.match(r"PAIR (\d+) ([\d.eE+-]+)", line)
        if m:
            times[int(m.group(1))] = float(m.group(2))
            cur = int(m.group(1))
        elif line.startswith("RESULT "):
            results[cur] = line[len("RESULT "):].strip()
    return times, results


def parse_maple(expr_str, n):
    """Return {perm tuple: sympy coefficient in T_i} from a Maple linear combination of X[...]."""
    Ts = {f"T[{i}]": sympy.Symbol(f"T{i}") for i in range(1, n)}
    s = expr_str
    for k, v in Ts.items():
        s = s.replace(k, str(v))
    perms = {}
    def repl(m):
        p = tuple(int(x) for x in m.group(1).split(","))
        name = "X_" + "_".join(map(str, p))
        perms[name] = p
        return name
    s = re.sub(r"X\[([\d, ]+)\]", repl, s)
    e = sympy.sympify(s)
    out = {}
    for name, p in perms.items():
        c = sympy.cancel(sympy.expand(e).coeff(sympy.Symbol(name)))
        if c != 0:
            out[p] = c
    return out


def schubmult_product(u, v, n, equivariant):
    R = DGx if equivariant else Gx
    prod = R(list(u)) * R(list(v))
    if equivariant:
        prod = prod.simplify(factor=False)
    out = {}
    for w, c in prod.items():
        w = tuple(w)
        if len(w) > n or max(w) > n:
            continue
        w = tuple(list(w) + list(range(len(w) + 1, n + 1)))
        c = sympy.sympify(str(c)).subs(sympy.Symbol("β"), -1) if "β" in str(c) else sympy.sympify(str(c))
        c = c.subs(sympy.Symbol("beta"), -1)
        out[w] = c
    return out


def to_T_form(coeff_y, n):
    """Substitute y_i -> 1 - t_i with T_i = t_i/t_{i+1}, i.e. compare in variables t_i = 1 - y_i."""
    subs = {sympy.Symbol(f"y_{i}"): 1 - sympy.Symbol(f"t{i}") for i in range(1, n + 1)}
    return sympy.cancel(coeff_y.subs(subs))


def maple_to_t(coeff_T, n):
    subs = {sympy.Symbol(f"T{i}"): sympy.Symbol(f"t{i}") / sympy.Symbol(f"t{i+1}") for i in range(1, n)}
    return sympy.cancel(coeff_T.subs(subs))


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    npairs = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    equivariant = "--equivariant" in sys.argv
    random.seed(1)
    maxlen = n * (n - 1) // 2
    def sample():
        while True:
            p = tuple(random.sample(range(1, n + 1), n))
            if 2 <= Permutation(list(p)).inv <= maxlen // 2:
                return p
    pairs = [(sample(), sample()) for _ in range(npairs)]

    t0 = time.perf_counter()
    mtimes, mres = run_maple(maple_script(n, pairs, equivariant))
    maple_wall = time.perf_counter() - t0

    stimes = []
    mismatches = 0
    for k, (u, v) in enumerate(pairs):
        t1 = time.perf_counter()
        sres = schubmult_product(u, v, n, equivariant)
        stimes.append(time.perf_counter() - t1)
        mparsed = parse_maple(mres[k], n)
        ok = set(mparsed) == set(sres)
        if ok:
            for w in sres:
                lhs = to_T_form(sres[w], n)
                rhs = maple_to_t(mparsed[w], n)
                if sympy.simplify(lhs - rhs) != 0:
                    ok = False
                    break
        if not ok:
            mismatches += 1
            print(f"MISMATCH {u} * {v}:\n  maple:     {mres[k]}\n  schubmult: {sres}")
        print(f"pair {k}: {u} * {v}  terms={len(sres)}  maple {mtimes[k]:.3f}s  schubmult {stimes[-1]:.3f}s  {'OK' if ok else 'FAIL'}")
    print(f"\nn={n} {'K_T' if equivariant else 'K'}: {len(pairs)} products, {mismatches} mismatches")
    print(f"maple total (compute only) {sum(mtimes.values()):.2f}s, wall incl. startup {maple_wall:.1f}s; schubmult total {sum(stimes):.2f}s")


if __name__ == "__main__":
    main()
