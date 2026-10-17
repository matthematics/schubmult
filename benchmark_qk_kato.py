"""Kato's projection QK(Fl_n) -> QK(Gr(m,n)) vs EquivCalc comin_qkmult / comin_qktmult.

Kato (arXiv:1906.09343): O^w -> O^{w W_P}, q_i -> 1 (i != m), q_m -> q is a ring homomorphism.
schubmult side: grothmult_q / grothmult_q_double, projected by ``apply_kato`` (beta-homogeneous
form), then beta = -1.
Usage: python benchmark_qk_kato.py m n [--equivariant]
"""
import itertools
import re
import subprocess
import sys
import time

import sympy

from schubmult import GeneratingSet, Permutation
from schubmult.mult.groth_double import normalize_coeff
from schubmult.mult.groth_quantum import grothmult_q
from schubmult.mult.groth_quantum_double import apply_kato, grothmult_q_double

MAPLE = "/mnt/c/Program Files/Maple 2026/bin.X86_64_WINDOWS/cmaple.exe"
EQUIVCALC = "C:/Users/matth/equivcalc/EquivCalc-1.0.1/equivcalc"
SCRIPT = "/mnt/c/Users/matth/equivcalc/qk_bench_tmp.mpl"


def redword(perm):
    w = list(perm)
    word = []
    while True:
        for i in range(len(w) - 1):
            if w[i] > w[i + 1]:
                w[i], w[i + 1] = w[i + 1], w[i]
                word.append(i + 1)
                break
        else:
            break
    return word[::-1]


def grassmannian_perms(m, n):
    out = []
    for S in itertools.combinations(range(1, n + 1), m):
        rest = [i for i in range(1, n + 1) if i not in S]
        out.append(tuple(list(S) + rest))
    return out


def run_maple(m, n, pairs, equivariant):
    fn = "comin_qktmult" if equivariant else "comin_qkmult"
    lines = [f'read "{EQUIVCALC}":', "with(equivcalc):", f"Gr({m},{n}):"]
    for k, (u, v) in enumerate(pairs):
        lines += [
            f"u := redexp_weyl({redword(u)}):",
            f"v := redexp_weyl({redword(v)}):",
            f"st := time(): r := {fn}(u, v): et := time() - st:",
            f'printf("PAIR {k} %a\\n", et):',
            'printf("RESULT %a\\n", weyl_sp(r)):',
        ]
    lines.append("quit;")
    with open(SCRIPT, "w") as f:
        f.write("\n".join(lines) + "\n")
    with open(SCRIPT) as f:
        out = subprocess.run([MAPLE, "-q"], stdin=f, capture_output=True, text=True, timeout=36000).stdout
    times, results = {}, {}
    for line in out.splitlines():
        mm = re.match(r"PAIR (\d+) ([\d.eE+-]+)", line)
        if mm:
            cur = int(mm.group(1))
            times[cur] = float(mm.group(2))
        elif line.startswith("RESULT "):
            results[cur] = line[len("RESULT "):].strip()
    return times, results, out


def parse_maple(expr_str, n):
    s = expr_str
    for i in range(1, n):
        s = s.replace(f"T[{i}]", f"T{i}")
    perms = {}

    def repl(mm):
        p = tuple(int(x) for x in mm.group(1).split(","))
        # Gr mode displays w_0 w w_0 relative to Fl mode (where X[w] <-> G_w)
        N = len(p)
        p = tuple(N + 1 - p[N - i] for i in range(1, N + 1))
        name = "X_" + "_".join(map(str, p))
        perms[name] = p
        return name

    s = re.sub(r"X\[([\d, ]+)\]", repl, s)
    e = sympy.expand(sympy.sympify(s))
    out = {}
    for name, p in perms.items():
        c = sympy.cancel(e.coeff(sympy.Symbol(name)))
        if c != 0:
            out[p] = c
    return out


def maple_to_t(coeff_T, n):
    subs = {sympy.Symbol(f"T{i}"): sympy.Symbol(f"t{i}") / sympy.Symbol(f"t{i+1}") for i in range(1, n)}
    return sympy.cancel(coeff_T.subs(subs))


def schubmult_to_t(coeff_y, n):
    subs = {sympy.Symbol(f"y_{i}"): 1 - sympy.Symbol(f"t{i}") for i in range(1, n + 1)}
    return sympy.cancel(coeff_y.subs(subs))


def schubmult_kato(u, v, m, n, equivariant):
    q = sympy.Symbol("q")
    if equivariant:
        y = GeneratingSet("y")
        d = grothmult_q_double({Permutation(list(u)): 1}, Permutation(list(v)), y, y)
    else:
        d = grothmult_q({Permutation(list(u)): 1}, Permutation(list(v)))
    d = apply_kato(d, [i for i in range(1, n) if i != m], n=n)
    if equivariant:
        d = {w: normalize_coeff(c, y) for w, c in d.items()}
    out = {}
    for w, c in d.items():
        c = sympy.sympify(str(c)).subs(sympy.Symbol("β"), -1).subs(sympy.Symbol("beta"), -1).subs(sympy.Symbol("q_1"), q)
        w = tuple(w)
        w = w + tuple(range(len(w) + 1, n + 1))
        out[w] = c
    return {w: c for w, c in out.items() if c != 0}


def main():
    m, n = int(sys.argv[1]), int(sys.argv[2])
    equivariant = "--equivariant" in sys.argv
    W = grassmannian_perms(m, n)
    pairs = [(u, v) for u in W for v in W if u <= v]
    t0 = time.perf_counter()
    times, results, raw = run_maple(m, n, pairs, equivariant)
    print(f"maple: {len(results)} products, compute {sum(times.values()):.2f}s, wall {time.perf_counter()-t0:.2f}s")
    mism = 0
    tsm = 0.0
    for k, (u, v) in enumerate(pairs):
        if k not in results:
            print("no maple result for", u, v)
            mism += 1
            continue
        ref = parse_maple(results[k], n)
        t0 = time.perf_counter()
        ours = schubmult_kato(u, v, m, n, equivariant)
        tsm += time.perf_counter() - t0
        ok = set(ref) == set(ours)
        if ok:
            for w in ref:
                lhs = maple_to_t(ref[w], n) if equivariant else ref[w]
                rhs = schubmult_to_t(ours[w], n) if equivariant else ours[w]
                if sympy.simplify(lhs - rhs) != 0:
                    ok = False
                    break
        if not ok:
            mism += 1
            print(f"MISMATCH {u} * {v}")
            print("  maple:", ref)
            print("  ours :", ours)
    print(f"Gr({m},{n}) {'QK_T' if equivariant else 'QK'}: {len(pairs)} products, {mism} mismatches; schubmult {tsm:.2f}s")


if __name__ == "__main__":
    main()
