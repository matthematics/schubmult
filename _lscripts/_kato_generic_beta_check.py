"""Generic-beta sanity checks for ``apply_kato``: for all pairs of minimal coset representatives in
S_n for the given block sizes, the projected product must be polynomial in beta and homogeneous of
degree l(u) + l(v) with deg beta = -1, deg G_w = l(w), deg q_i = (size of block i) + (size of block i+1).
Usage: python _lscripts/_kato_generic_beta_check.py n1 n2 ... [--equivariant]
"""
import itertools
import sys

import sympy

from schubmult import GeneratingSet, Permutation
from schubmult.mult.groth_double import normalize_coeff
from schubmult.mult.groth_quantum import grothmult_q
from schubmult.mult.groth_quantum_double import apply_kato, grothmult_q_double


def main():
    blocks = [int(a) for a in sys.argv[1:] if not a.startswith("--")]
    equivariant = "--equivariant" in sys.argv
    n = sum(blocks)
    parabolic_index, start = [], 0
    for size in blocks:
        parabolic_index += list(range(start + 1, start + size))
        start += size
    qdeg = {sympy.Symbol(f"q_{k + 1}"): blocks[k] + blocks[k + 1] for k in range(len(blocks) - 1)}
    beta = sympy.Symbol("β")
    y = GeneratingSet("y")
    reps = [Permutation(list(p)) for p in itertools.permutations(range(1, n + 1)) if not any(d + 1 in parabolic_index for d in Permutation(list(p)).descents())]
    bad = 0
    total = 0
    for u, v in itertools.combinations_with_replacement(reps, 2):
        total += 1
        if equivariant:
            d = grothmult_q_double({u: 1}, v, y, y)
        else:
            d = grothmult_q({u: 1}, v)
        d = apply_kato(d, parabolic_index, n=n)
        target = u.inv + v.inv
        for w, c in d.items():
            c = sympy.sympify(str(normalize_coeff(c, y)))
            if c == 0:
                continue
            num, den = sympy.fraction(c)
            ok = True
            # denominators must be products of (1 + beta y_i): no pure beta power survives at y = 0
            if den.subs({s: 0 for s in den.free_symbols if str(s).startswith("y_")}) != 1:
                ok = False
            for term in sympy.Add.make_args(sympy.expand(num)):
                pd = term.as_powers_dict()
                deg = -int(pd.get(beta, 0)) + sum(int(pd.get(q, 0)) * dq for q, dq in qdeg.items())
                deg += sum(int(e) for s, e in pd.items() if str(s).startswith("y_"))
                if deg + Permutation(w).inv != target:
                    ok = False
            if den != 1:
                for term in sympy.Add.make_args(sympy.expand(den)):
                    pd = term.as_powers_dict()
                    deg = -int(pd.get(beta, 0)) + sum(int(e) for s, e in pd.items() if str(s).startswith("y_"))
                    if deg != 0:
                        ok = False
            if not ok:
                bad += 1
                print(f"BAD {u} * {v} -> {w}: {c}")
    print(f"blocks {blocks} {'QK_T' if equivariant else 'QK'}: {total} products, {bad} bad coefficients")


if __name__ == "__main__":
    main()
