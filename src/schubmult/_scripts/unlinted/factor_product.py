from schubmult import *
from sympy import pretty_print, S
from schubmult.utils.schub_lib import elem_sym_perms_op, elem_sym_perms
from schubmult.combinatorics.crystal_graph import CrystalGraphTensor
from schubmult.symbolic import expand_seq
from schubmult.rings.polynomial_algebra import *

def the_prod(tup):
    if len(tup) == 0:
        return RCGraph()
    return the_prod(tup[:-1]).resize(len(tup[-1])).squash_product(tup[-1])

if __name__ == "__main__":
    import sys
    import itertools

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    r = DualRCGraphRing()
    factorizations = {RCGraph(): set([((),())])}
    for k in range(1, n):
        new_factorizations = {}
        for p in range(0, k + 1):
            for rc, factor_set in factorizations.items():
                rc = rc.resize(k)
                for elem_rc in RCGraph.all_rc_graphs(uncode([0] * (k - p) + [1] * p), k):
                    for new_rc in r(rc) * r(elem_rc):
                        new_factorizations[new_rc] = new_factorizations.get(new_rc, set())
                        for (factor_tuple, weight) in factor_set:
                            new_factorizations[new_rc].add(((*factor_tuple, elem_rc),(*weight, elem_rc.perm.inv)))
                            break
        factorizations = new_factorizations

    print("TEST FK FACTORIZATON")

    
    for base_perm in perms:
        beepy = 0
        buildup = 1
        for rc_b in RCGraph.all_rc_graphs(base_perm, n - 1):
            for perm in perms:
                def word_elem(p, k, *args):  # noqa: ARG001
                    import numpy as np
                    numvars = n - 1
                    if p > k or p < 0:
                        return PA.zero
                    if p == 0:
                        return PA.one
                    vec = np.zeros(numvars, dtype=int)
                    vec[k - 1] = p
                    return PA.from_dict({tuple(vec.tolist()): 1})
                        
                
                new_buildup = 0
                for rc in RCGraph.all_rc_graphs(perm, n - 1):
                    
                    for (factor_tuple, weight) in factorizations[rc]:
                        for (factor_tuple2, weight2) in factorizations[rc_b]:
                            term = 1
                            # if weight != weight2:
                            #     continue
                        # sem = Sx(perm).in_SEM_basis(elem_func=word_elem)
                        # neg_coeff = sem.get(weight, 0)
                        # * FA(*weight).change_basis(SchubertBasis).get(perm * Permutation.w0(n), n - 1)
                        
                        # if neg_coeff == 0:
                        #     print(f"Pantalones park of potato {neg_coeff=} {rc.perm=} {weight=} {sem=}")
                        #     continue
                            #perm_dict = {(RCGraph(), Permutation([])): 1}
                            
                            vec_sum = ()
                            
                            
                            for k, (factor, factor2) in enumerate(zip(factor_tuple,factor_tuple2), start=1):
                                new_perm_dict = 0
                                #the_diff = sum(factor.length_vector)
                                
                                    # for up_perm, diff in elem_sym_perms(p, k, k):
                                    #     if diff != the_diff:
                                    #         continue
                                    #     # if len(up_perm) > k + 1:
                                    #     #      continue
                                    #     #test_vec = tuple([1 if up_perm[i] != p[i] else 0 for i in range(k)])
                                    #     #crapeff = ASx(up_perm, k).change_basis(WordBasis).get(vec_sum, 0)
                                    #     #if crapeff > 0:
                                term *=   (r(factor) * r(factor2)) @ (Sx(factor.perm)*Sx(factor2.perm))
                            new_buildup += buildup * term
                            sparse_termo = set(new_buildup.keys())
                            for (rcc, ppp) in sparse_termo:
                                if rcc.perm != ppp:
                                    del new_buildup[(rcc, ppp)]
                    buildup += new_buildup
                            #buildup += (r@Sx).from_dict(perm_dict)
                    #buildup += (r@Sx).from_dict(perm_dict)
                
            try:
                potato = Sx(base_perm) * Sx(perm)
                #buildup = (r@Sx).from_dict(perm_dict)
                #potato = Sx.from_dict({k: v for k, v in (Sx(base_perm) * Sx(perm)).items() if len(k) <= n})
                potato2 = sum([coeff * Sx(up_perm) for (rc, up_perm), coeff in buildup.items() if rc.is_principal and rc.perm == up_perm])
                assert potato2.almosteq(potato), f"Failed on {base_perm} * {perm}\n{buildup}\n{potato}"
                print(f"Passed for {base_perm} * {perm}")
            except AssertionError as e:
                print(e  )
                
            
    print("Potanto")
    sys.exit(0)
    for perm in perms:
        for rc in RCGraph.all_rc_graphs(perm, n - 1):
            factor_tuple = factorizations[rc]
            tensor = factor_tuple[0]
            for i in range(1, len(factor_tuple)):
                tensor = CrystalGraphTensor(tensor, factor_tuple[i].resize(n-1))
            tensor, _ = tensor.to_lowest_weight()
            rc_lw = rc.to_lowest_weight()[0]
            assert rc_lw.crystal_weight == tensor.crystal_weight, f"Failed on {rc}\n{tensor}\n{rc_lw}"
            to_tuple = []
            the_factor = tensor
            for i in range(len(factor_tuple) - 1):
                to_tuple.insert(0, the_factor.factors[1])
                the_factor = the_factor.factors[0]
            to_tuple.insert(0, the_factor)
            product = RCGraph()
            for i, rct in enumerate(to_tuple):
                rct = rct.resize(i + 1)
                product = product.resize(len(rct)).squash_product(rct)
            assert product == rc.to_lowest_weight()[0], f"Failed on {rc}\n{to_tuple}\n{product}"

            tensor, _ = tensor.to_highest_weight()

            to_tuple = []
            the_factor = tensor
            for i in range(len(factor_tuple) - 1):
                to_tuple.insert(0, the_factor.factors[1])
                the_factor = the_factor.factors[0]
            to_tuple.insert(0, the_factor)
            product = RCGraph()
            for i, rct in enumerate(to_tuple):
                rct = rct.resize(i + 1)
                product = product.resize(len(rct)).squash_product(rct)
            assert product == rc.to_highest_weight()[0], f"Failed on {rc}\n{to_tuple}\n{product}"
            # assert product.to_highest_weight()[0] == rc.to_highest_weight()[0], f"Failed on {rc}\n{to_tuple}\n{product}"
                
            