from schubmult import *
from schubmult.rings.combinatorial.dual_rc_graph_ring import DualRCGraphRing
from schubmult.rings.polynomial_algebra import *
from sympy import pretty_print
import argparse

r = DualRCGraphRing()

def word_elem(p, k, *args):  # noqa: ARG001
    import numpy as np
    #numvars = n - 1
    if p > k or p < 0:
        return PA.zero
    if p == 0:
        return 1
    vec = np.zeros(10, dtype=int)
    vec[k - 1] = p
    return PA.from_dict({tuple(vec.tolist()): 1})

                    

def dual_elem_sym(p, k, nn, varvar):
    if varvar != 1:
        ret = Sx.zero @ r.zero
        for i in range(p + 1):
            ret += Sx.from_expr(varvar**i) @ r.schub((uncode([0] * (k - p + i)+[1] * (p - i))), nn)
    else:
        return r.schub((uncode([0] * (k - p)+[1] * (p))), nn)
        #Sx(uncode([0] * (k - p)+[i])) 
    #pretty_print(result)
    return ret

def dual_elem_sym2(p, k, nn, varvar):
    ret = FA() @ r.zero
    for i in range(p + 1):
        ret += FA(i) @ r.schub((uncode([0] * (k - p + i)+[1] * (p - i))), nn)
        #Sx(uncode([0] * (k - p)+[i])) 
    #pretty_print(result)
    return ret

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute terms for squash_cauchy with optional custom product order.",
    )
    parser.add_argument("n", type=int, help="Positive integer n.")
    parser.add_argument(
        "--order",
        nargs="+",
        type=int,
        help="Space-separated order of indices 1..n-1 indicating product order.",
    )
    args = parser.parse_args()

    n = args.n
    if n < 2:
        parser.error("n must be at least 2")

    default_order = list(range(1, n))
    if args.order is None:
        product_order = default_order
    else:
        product_order = args.order
        expected = set(default_order)
        actual = set(product_order)
        if len(product_order) != len(default_order) or actual != expected:
            parser.error(f"--order must contain each integer from 1 to {n - 1} exactly once")
    
    perms = Permutation.all_permutations(n)
    for perm in perms:
        if perm.inv == 0:
            continue
        sem = Sx(perm).in_SEM_basis(elem_func=word_elem)
        # print(f"{sem=}")
        the_fart_potato = 0
        for weightburger, coeff in sem.items():
            if coeff == 0:
                continue
            term = None
            for k in product_order:
                a = weightburger[k - 1] if k - 1 < len(weightburger) else 0
                factor = dual_elem_sym(a, k, n, 1)
                if term is None:
                    term = coeff * factor
                else:
                    term = term @ factor
            the_fart_potato += term
        pretty_print(the_fart_potato)
        #assert all(rc.perm == perm for rc in the_fart_potato), f"Test failed for {perm}, got {the_fart_potato}"
            #print(f"{perm=}, {weightburger=}, {coeff=}")
    raise SystemExit(0)
    result0 = 1
    for i in range(1, n):
        result0 *= dual_elem_sym2(i, i, n - 1, Sx.genset[n - i])

    result = None
    for i in range(1, 2 * (n - 1) + 1):
        j = (i + 1) //2
        if result is None:
            result = dual_elem_sym(j, j, n, Sx.genset[2*(n - 1) + 1 - i])
        else:
            result = result * dual_elem_sym(j, j, n, Sx.genset[2*(n - 1) + 1 - i])
    dct_result0_spank = {}
    dct_result0 = {}
    dct_result0_schub = {}
    for (u, v), coeff in result0.items():
        dct_result0[u] = dct_result0.get(u, 0) + coeff * r(v)

    for (u, v), coeff in result.items():
        dct_result0_spank[u] = dct_result0_spank.get(u, 0) + coeff * r(v)

    result0_schub = 0
    for (u, v), coeff in result0.items():
        result0_schub += coeff * FA(*u).change_basis(SchubertBasis) @ r(v)
        #dct_result0[u].add(v)
    
    for (u, v), coeff in result0_schub.items():
        dct_result0_schub[u]    = dct_result0_schub.get(u, 0) + coeff * r(v)
    # pretty_print(result)
    # # square
    # # new_result = 0
    # # for (u, v), coeff in result.items():
    # #     for (u2, v2), coeff2 in result.items():
    # #         new_result += coeff * coeff2 * Sx(u) *Sx(u2) @ r(v) @ r(v2)
    # # pretty_print(new_result)
    # # coprod n - 1

    # for (u, v), coeff in result.items():
    #     coprod = Sx(u).coproduct(*list(range(1,n//2 + 1)))
    #     for (u1, u2), coeff2 in coprod.items():
    #         #print(f"{u1=}, {u2=}, {v=}, {coeff=}, {coeff2=}")
    #         pretty_print(coeff * coeff2 * Sx(u1) @ Sx(u2) @ r(v))

    def cope(rc_graph, result):
        ret = 0
        for (u, rc), coeff in result.items():
            if rc != rc_graph:
                continue
            coprod = Sx(u).coproduct(*list(range(1, 2*(n - 1) + 1, 2)))
            for (u1, u2), coeff2 in coprod.items():
                u1 = Permutation(u1)
                u2 = Permutation(u2)
                if u1 not in dct_result0_spank or u2 not in dct_result0_spank:
                    continue
                for rc1, coeff_rc1 in dct_result0_spank[u1].items():
                    for rc2, coeff_rc2 in dct_result0_spank[u2].items():
                        #if tuple([a + b for a, b in zip(rc1.length_vector, rc2.length_vector)]) == rc_graph.length_vector:
                        ret += coeff * coeff2 *  r(rc1) @ ASx(rc2.perm, len(rc2))
                            #ret += coeff * coeff2 * r(rc1) @ FA(*rc2.length_vector)
                            #ret += coeff * coeff2 * r(rc1) @ ASx(rc2.perm, len(rc2))
                            #break
        return ret
                #ret += coeff * coeff2 * ASx(u1 * w0, n - 1) @ ASx(u2, n - 1) @ r(rc)
                #ret += coeff * coeff2 * Sx(u1 * w0) @ Sx(u2) @ r(rc)
    
    def cope2(rc_graph, dct_result0):
        ret = 0
        for (word, rc) in result0:
            if rc != rc_graph:
                continue
            cop = FA(*word).coproduct()
            for (word1, word2), coeff in cop.items():
                for rc1, coeff1 in dct_result0[word1].items():
                    for rc2, coeff2 in dct_result0[word2].items():
                        if tuple([a + b for a, b in zip(rc1.length_vector, rc2.length_vector)]) == rc_graph.length_vector:
                            ret += coeff * coeff1 * coeff2 * r(rc1) @ FA(*word2)
        return ret

    def cope2_asx(rc_graph):
        ret = 0
        for (word, rc) in result0_schub:
            if rc != rc_graph:
                continue
            cop = ASx(*word).coproduct()
            for (word1, word2), coeff in cop.items():
                for rc1, coeff1 in dct_result0_schub[word1].items():
                    for rc2, coeff2 in dct_result0_schub[word2].items():
                        if tuple([a + b for a, b in zip(rc1.length_vector, rc2.length_vector)]) == rc_graph.length_vector:
                            ret += coeff * coeff1 * coeff2 * r(rc1) @ ASx(*word2)
        return ret

    perms = Permutation.all_permutations(n)
    for perm in perms:
        pants = 0   
        pofer = 0
        for rc in RCGraph.all_rc_graphs(perm, n - 1):
            #result = r(rc)
            # pretty_print(rc)
            
                #pretty_print(cop)
            #cop = cope2(rc, dct_result0)
            cop = cope2_asx(rc)
            if cop == 0:
                continue
            pofer += 1
            # test comodule condition
            comod_elem1 = 0
            for (rc1, word1), coeff1 in cop.items():
                #comod_elem1 += coeff1 * r(rc1)@ (FA(*word1).coproduct())
                comod_elem1 += coeff1 * r(rc1)@ (ASx(*word1).coproduct())

            comod_elem2 = 0
            for (rc1, word1), coeff1 in cop.items():
                # comod_elem2 += coeff1 * cope2(rc1, dct_result0) @ FA(*word1)
                comod_elem2 += coeff1 * cope(rc1, result) @ ASx(*word1)
            try:
                assert comod_elem1.almosteq(comod_elem2), f"Comodule condition failed for {rc} with coproduct {cop}\n{comod_elem1-comod_elem2}"
            except AssertionError as e:
                #print(e)
                # pretty_print(comod_elem1)
                # pretty_print(comod_elem2)
                pants += 1
        print(f"Checked comodule condition for all RC graphs of {perm}, failed {pants} times out of {pofer}")
        #print("-----")