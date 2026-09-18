from schubmult import *
from schubmult.abc import *
from schubmult.rings.polynomial_algebra import PA
Sy = SingleSchubertRing(y)

def dom_from_seq(seq2):
    N = len(seq2)
    cd = [0] * N
    for i in range(1, len(seq2) + 1):
        for j in range(i):
            cd[j] += seq2[i - 1]
    return uncode(cd)

def cauchy_tensor(perm, n):
    the_dub = DSx(perm).expand()
    first_result = Sx.from_expr(the_dub)
    result = (Sx@PA).zero
    for perm, coeff in first_result.items():
        schub = PA.from_expr(coeff, length = n - 1)#Sy.from_expr(coeff)
        result += Sx(perm) @ schub
    return result

if __name__ == "__main__":
    import sys
    from sympy import pretty_print
    from schubmult.utils.perm_utils import artin_sequences
    from schubmult.abc import *
    from schubmult.rings.polynomial_algebra import PA
    from schubmult.rings.schubert import SingleSchubertRing

    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)
    seqs = artin_sequences(n)
    for perm in perms:
        #for seq in seqs:
        seq = tuple([1] * (n - 1))
            # if sum(seq) != perm.inv:
            #     continue
            #rev_seq = tuple(reversed(seq))
            # cop = FA(*seq).coproduct()
            # # treat as dominant
        result = 0
        copr = FA(*seq).coproduct()
        dom = dom_from_seq(seq)
        tensor0 = cauchy_tensor(dom, n)
        perminv = (perm) * dom
        for (word1, word2), coeff in copr.items():
            dom1 = dom_from_seq(word1, n)
            dom2 = dom_from_seq(word2, n)
            tensor1 = cauchy_tensor(dom1, n)
            tensor2 = cauchy_tensor(dom2, n)
            for (perm11, bangtash1), coeff1 in tensor1.items():
                for (perm21, bangtash2), coeff2 in tensor2.items():
                    if perm12 == perminv and perm22 == perminv:
                        result += coeff * coeff1 * coeff2 *Sx(perm11) @ Sx(perm21)
            print(result)
            #print(cauchy_tensor(dom))
            # schubs = Sx.from_expr(DSx(~dom).expand().subs({y[i]: 1 for i in range(len(y))}))
            # for schub in schubs:
            # part = tuple(dom.trimcode)
            # cems = RCGraph.full_CEM(schub * dom, n, part)
            # for rc, cem_dict in cems.items():
            #     for tup, coeff in cem_dict.items():
            #         if coeff != 0:
            #             result += coeff * FA(*tup) @ CompositionSchubertBasis
            # # result = 0
            # # def flip_sign_vars(expr, vars):
            # #     subs_dict = {vars[i]: -vars[i] for i in range(len(vars))}
            # #     return expr.subs(subs_dict)
            
            # # def sx_potato(elem, ring2, word):
            # #     bong = Sx.from_expr(elem.expand())
            # #     paper = (FA @ ring2).zero
            # #     for perm, coeff in bong.items():
            # #         coeff = flip_sign_vars(coeff, ring2.genset)
            # #         stork = PA.from_expr(Sx(perm).expand(), length=len(word))
            # #         boing = stork.get(tuple(word), 0)
            # #         if boing != 0:
            # #             paper += boing * FA(*word) @ ring2.from_expr(coeff)
                    
            # #     return paper

            # # yring = SingleSchubertRing(y)
            # # zring = SingleSchubertRing(z)
            # # for (word1, word2), coeff in cop.items():
            # #     prd1 = 1
            # #     #prd2 = 1
            # #     for index, degree in enumerate(word1, start=1):
            #         prd1 *= Sx.from_expr(e(degree, n + 1 - index, x[1:]))
            #     # for index, degree in enumerate(word2, start=1):
            #     #     prd2 *= Sx.from_expr(e(n + 1 - index - degree, n + 1 - index, x[1:]))
            #     #dom1 = dom_from_seq(word1)
            #     #dom2 = dom_from_seq(word2)
            #     #result += coeff * sx_potato(DSx(dom1).expand(), yring, word1) @ sx_potato(DSx(dom2, "z").expand(), zring, word2)
            #     #result += Sx.from_expr(e(n + 1 - ))
            #     result += coeff * prd1 @ FA(*word2).change_basis(CompositionSchubertBasis)
            # print(f"Seq: {seq}, Result: {result}")
