from schubmult import *
from schubmult.combinatorics.double_forest import long_word, sylvester_word, forest_from_code, canonical_forest_from_word, _max_leaf, letter_weight
from schubmult.combinatorics.indexed_forests import letterpair, omega_insertion
from schubmult.rings.polynomial_algebra import *
from schubmult.symbolic import expand

def vine_subword_polynomial(code, x_gen, t_gen):
    """Vine subword model from arXiv:2504.15234 (Bergeron-Gagnon-Nadeau-Spink-Tewari).

    Iterates subwords pi of the long word `long_word(n)` of size
        sz = sum(code) = ell(uncode(code))
    and sums wt(pi) under one of two filters on the *value sequence*
    (treating barred and unbarred letters by their value only):

      mode='forest'   :  values must be Sylvester-equivalent to sylvester_word(F),
                         i.e. lie in Syl(F).  Equals P_F(x;t)  (Theorem 5.1).
      mode='schubert' :  values must be a reduced word for perm=uncode(code),
                         i.e. lie in Red(perm).  Equals S_perm(x;t) (Theorem 6.1).

    Weights (Sylvester column convention, paper p.19):
        wt(j^(k))      = x_k - t_j
        wt(barred j^(k)) = t_j - t_k
    realised by `letter_weight` with bracket-indexed `x_gen`, `t_gen`.
    """
    from schubmult import Permutation, uncode
    from itertools import combinations

    code = list(code)
    perm = uncode(code)
    sz = sum(code)
    assert sz == perm.inv, f"sum(code)={sz} != perm.inv={perm.inv}"

    F = forest_from_code(code)
    max_leaf = 0
    for tt in F:
        max_leaf = max(max_leaf, _max_leaf(tt))
    n = max(max_leaf - 1, sz, len(perm.trimcode), 1)
    
    word = long_word(n)

    def subwords():
        stack = [(perm, [], word)]
        while len(stack) > 0:
            working_perm, working_subword, remaining_word = stack.pop()
            if working_perm.inv == 0:
                yield working_subword
                continue
            if len(remaining_word) < working_perm.inv:
                continue
            descent = remaining_word[-1]["value"]
            new_remaining_word = remaining_word[:-1]
            if working_perm[descent - 1] > working_perm[descent]:
                new_perm = working_perm.swap(descent - 1, descent)
                stack.append((new_perm,  [remaining_word[-1]] + working_subword, new_remaining_word))
            stack.append((working_perm, working_subword, new_remaining_word))

    def _pair_label(values):
        counts = {}
        out = []
        for a in values:
            counts[a] = counts.get(a, 0) + 1
            out.append(letterpair(a, counts[a]))
        return tuple(out)

    sw = sylvester_word(F)
    sylv_target_canon = canonical_forest_from_word(sw)
    sylvester_letters = None
    for idx in combinations(range(len(word)), sz):
        vs = [word[i]["value"] for i in idx]
        if canonical_forest_from_word(vs) == sylv_target_canon:
            sylvester_letters = tuple(word[i] for i in idx)
            break

    total = 0
    for subword in subwords():
        values = tuple(letter["value"] for letter in subword)
        w = 1
        lbs, dec = omega_insertion(_pair_label(tuple(reversed(values))))
        for nd in dec.inorder_traversal:
            w = w * letter_weight(sylvester_letters[nd[1] - 1], x_gen, t_gen)
        total = total + w
    return expand(total)


if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    comps = [perm.pad_code(n - 1) for perm in perms]
    x, t = GeneratingSet("x"), GeneratingSet("t")
    DForest = PolynomialAlgebra(DoubleForestPolyBasis(x,t))

    for comp in comps:
        subval = vine_subword_polynomial(comp, x, t)
        forval = DForest(*comp).expand()
        diff = (subval - forval).expand()
        assert diff==0, f"Mismatch for comp={comp}\n{diff=}\n{subval=}\n{forval=}"
        print(f"Verified vine_subword_polynomial for comp={comp}\n{diff=}")
    # import argparse

    # parser = argparse.ArgumentParser(description="Compute vine subword polynomial")
    # parser.add_argument("code", type=str, help="Code as a comma-separated list of integers")
    # parser.add_argument("--n", type=int, default=None, help="Optional n value")
    # parser.add_argument("--mode", type=str, choices=["forest", "schubert"], default="forest", help="Mode: 'forest' or 'schubert'")
    # args = parser.parse_args()

    # code = [int(x) for x in args.code.split(",")]
    # result = vine_subword_polynomial(code, x_gen=var("x"), t_gen=var("t"), n=args.n, mode=args.mode)
    # print(result)