import math
import sys

from schubmult import ASx, Permutation, Sx, uncode
from schubmult.abc import x
from schubmult.rings import FA, WordBasis
from schubmult.rings.free_algebra_basis import SchubertBasis
from schubmult.rings.rc_graph_module import RCGraph, RCGraphModule, RCGraphTensor, TensorModule, all_fa_degree
from schubmult.symbolic import S, expand_seq
from schubmult.utils.perm_utils import artin_sequences


def product_too_big(perm1, perm2, n):
    return any(len(perm) > n for perm in (Sx(perm1)*Sx(perm2)).keys())

def vector_sum(v1, v2):
    return tuple([a + b for a, b in zip(v1, v2)])

def main():
    n = int(sys.argv[1])

    seqs = set()
    degree = (n * (n - 1)) // 2
    for deg in range(degree + 1):
        seqs.update(all_fa_degree(deg, n-1))

    unit_rc_module = RCGraphModule({RCGraph(): 1})

    # 100% positive!

    unit_tensor_rc_module = TensorModule.ext_multiply(unit_rc_module, unit_rc_module)

    solution_module = TensorModule()
    solution_module2 = TensorModule()
    solution_module3 = TensorModule()

    aseqs = artin_sequences(n - 1)

    perms = Permutation.all_permutations(n)

    rc_graphs = {perm: RCGraph.all_rc_graphs(perm, n - 1) for perm in perms}

    rc_graphs_by_weight = {}

    for perm, rcs in rc_graphs.items():
        for rc in rcs:
            dct = rc_graphs_by_weight.get(perm, {})
            st = dct.get(rc.length_vector(),set())
            st.add(rc)
            dct[rc.length_vector()] = st
            rc_graphs_by_weight[perm] = dct
            

    principal_rcs = {perm: next(iter([rc for rc in rc_graphs[perm] if rc.is_principal])) for perm in perms}
    check3 = TensorModule()

    test_val = {}
    # for perm0 in perms:
    #     perm_elem = ASx(perm0, n-1).change_basis(WordBasis)
    #     for seq, coeff0 in perm_elem.items():
    #         if coeff0 == 0:
    #             continue
    #         rc_module = RCGraphModule()
    #         right_factor = FA(*seq).change_basis(SchubertBasis).coproduct()
    #         for weight in rc_graphs_by_weight[perm0]:
    #             rc_module += len(rc_graphs_by_weight[perm0][weight]) * coeff0 * RCGraphModule(dict.fromkeys(rc_graphs_by_weight[perm0][weight],1))
    #             # test_val[(perm0, weight)] = test_val.get((perm0, weight), 0) + len(rc_graphs_by_weight[perm0][weight]) * coeff0
    #         print(rc_module)
    #         solution_module3 += TensorModule.ext_multiply(rc_module, right_factor)

    # SAME FOR aLL RC GRAPHS BY WEIGHT
    # POSITIVE FOR PRINCIPAL
    # THEREFORE ALL POSITIVE
    for perm0 in perms:
        perm_elem = ASx(perm0, n-1).change_basis(WordBasis)
        for perm1 in perms:
            for perm2 in perms:
                for weight in rc_graphs_by_weight[perm0]:
                    the_sum = 0
                    weight_coeff = len(rc_graphs_by_weight[perm0][weight])
                    flag = False
                    for weight1 in rc_graphs_by_weight[perm1]:
                        for weight2 in rc_graphs_by_weight[perm2]:
                            if vector_sum(weight1, weight2) == weight:
                                flag = True
                            coeff0 = perm_elem.get(vector_sum(weight1, weight2), 0)
                            the_sum += weight_coeff * coeff0
                            for rc1 in rc_graphs_by_weight[perm1][weight1]:
                                for rc2 in rc_graphs_by_weight[perm2][weight2]:
                                    assert weight_coeff == FA(*weight).change_basis(SchubertBasis).get((perm0, n-1))
                                    rc_module = RCGraphModule(dict.fromkeys(rc_graphs_by_weight[perm0][weight],1))
                                    solution_module3 += weight_coeff * coeff0 *TensorModule.ext_multiply(rc_module, (ASx@ASx)(((perm1, n-1),(perm2, n-1))))
                    assert the_sum == 0 or (the_sum == 1 and flag)

    for (rc, ((perm1, _), (perm2, _))), coeff in solution_module3.items():
        assert coeff >= 0

    # print(test_val)

    # for (perm0, seq), coeff in test_val.items():
    #     assert coeff == Sx(expand_seq(seq, x)).get(perm0,0), f"Fail on {perm0} {seq} {coeff} {Sx(expand_seq(seq, x)).get(perm0,0)}"
    # exit()
    # for seq in aseqs:
    #     schub_elem = FA(*seq).change_basis(SchubertBasis)
    #     for perm0, coeff0 in schub_elem.items():
    #         loop_elem = ASx(*perm0) * unit_rc_module
    #         solution_module3 += coeff0 * TensorModule.ext_multiply(loop_elem,FA(*seq).coproduct() * unit_tensor_rc_module)

        #         for rc2 in rc_module:
        #             check3[rc] = check3.get(rc, TensorModule()) + coeff0 * (FA(*seq).coproduct() * unit_tensor_rc_module)

        #elem1 = FA(*seq).change_basis(SchubertBasis)
        #elem2 = elem1.coproduct()

        # for key, coeff1 in elem1.items():
        #     mod1 = ASx(*key) * unit_rc_module
        #     mod1 = RCGraphModule({rc: v for rc, v in mod1.items() if len(rc.perm) <= n})
        #     # for (k1, k2), coeff2 in elem2.items():
        #     #     mod2 = ASx(*k1) * unit_rc_module
        #     #     mod2 = RCGraphModule({rc: v for rc, v in mod2.items() if len(rc.perm) <= n})
        #     #     mod3 = ASx(*k2) * unit_rc_module
        #     #     mod3 = RCGraphModule({rc: v for rc, v in mod3.items() if len(rc.perm) <= n})
        #     solution_module2 += coeff1 * TensorModule.ext_multiply(mod1, FA(*seq).change_basis(SchubertBasis).coproduct() * unit_tensor_rc_module)

    # print(solution_module)
    # exit()

    # for rc, coeff in check3.items():
    #     print(f"{rc.perm.trimcode}")
    #     print(f"{coeff}")

    # exit()
    # THIS IS THE CORRECT COPRODUCT
    # for seq in seqs:
    #     perm = uncode(seq)
    #     elem = ASx(perm, n - 1).change_basis(WordBasis)
        
    #     for word_double, coeff_double in elem.items():
            


    #         mod = FA(*word_double) * unit_rc_module
    #         for rc, coeff2 in mod.items():
    #             assert rc.length_vector() == word_double
    #             for word, coeff in elem.items():
    #                 modmod = FA(*word).coproduct() * unit_tensor_rc_module
    #                 solution_module += coeff * coeff_double * coeff2 * TensorModule.ext_multiply(1 * rc, modmod)
    #                 solution_module2 += coeff * coeff_double * coeff2 * TensorModule.ext_multiply(RCGraphModule(dict.fromkeys(RCGraph.all_rc_graphs(rc.perm, n - 1),1)), modmod)

    coprods_interim = {}
    coprods = {}
    coprods_length_vector = {}
    coprods_length_vector2 = {}
    coprods_interim_rc = {}
    products = {}
    products_rc = {}
    products_by_weight = {}
    weight_sum = {}
    rc_sum = {}
    products_rc2 = {}
    # for (rc, (rc1, rc2)), coeff in solution_module.items():
    #     perm = rc.perm
    #     perm1, perm2 = rc1.perm, rc2.perm
    #     if len(perm1) > n or len(perm2) > n or len(perm) > n:
    #         continue
    #     # assert coeff >= 0 NOPE
    #     # coprods_interim[rc.perm] = coprods_interim.get(rc.perm, TensorModule()) + coeff * TensorModule.ext_multiply(1*rc1, ASx(perm2, len(rc2)))
    #     # coprods_interim_rc[rc.perm] = coprods_interim.get(rc.perm, TensorModule()) + coeff * TensorModule.ext_multiply(Sx(rc1.perm), 1 * rc2)
    #     # coprods_length_vector[rc.length_vector()] = coprods_length_vector.get(rc.length_vector(), 0) + coeff * (FA@FA)((rc1.length_vector(), rc2.length_vector()))
    #     # coprods_length_vector2[rc.length_vector()] = coprods_length_vector2.get(rc.length_vector(), 0) + coeff * (ASx@ASx)(((rc1.perm,len(rc1)), (rc2.perm,len(rc2))))
    #     products[(rc1.perm, rc2.perm)] = products.get((rc1.perm, rc2.perm), RCGraphModule()) + coeff * RCGraphModule(dict.fromkeys(RCGraph.all_rc_graphs(rc.perm, n - 1), 1))
        
    #     #products_by_weight[weight] = products_by_weight.get(weight, RCGraphModule()) + coeff * , 1))
    
    # for (rc, (rc1, rc2)), coeff in solution_module2.items():
    #     perm = rc.perm
    #     perm1, perm2 = rc1.perm, rc2.perm
    #     if len(perm1) > n or len(perm2) > n or len(perm) > n:
    #         continue
    #     # assert coeff >= 0 NOPE
    #     # coprods_interim[rc.perm] = coprods_interim.get(rc.perm, TensorModule()) + coeff * TensorModule.ext_multiply(1*rc1, ASx(perm2, len(rc2)))
    #     # coprods_interim_rc[rc.perm] = coprods_interim.get(rc.perm, TensorModule()) + coeff * TensorModule.ext_multiply(Sx(rc1.perm), 1 * rc2)
    #     # coprods_length_vector[rc.length_vector()] = coprods_length_vector.get(rc.length_vector(), 0) + coeff * (FA@FA)((rc1.length_vector(), rc2.length_vector()))
    #     # coprods_length_vector2[rc.length_vector()] = coprods_length_vector2.get(rc.length_vector(), 0) + coeff * (ASx@ASx)(((rc1.perm,len(rc1)), (rc2.perm,len(rc2))))
    #     products_rc[(rc1.perm, rc2.perm)] = products_rc.get((rc1.perm, rc2.perm), RCGraphModule()) + coeff * rc
    #     #products_by_weight[(perm1, perm2)] = products_by_weight.get((perm1, perm2), RCGraphModule()) + coeff * rc

    # products_by_bacon = {}
    # coproducts_by_bacon = {}
    # products_by_bacon2 = {}

    # for (rc, (rc1, rc2)), coeff in solution_module2.items():
    #     if len(rc.perm) > n or len(rc1.perm) > n or len(rc2.perm) > n:
    #         continue
    #     #perm1, perm2 = key1[0], key2[0]
    #     # perm = rc.perm
    #     # products_by_bacon[(rc1, rc2)] = products_by_bacon.get((rc1, rc2), 0) + coeff * Sx(rc.perm)
    #     #products_by_bacon2[(rc1, rc2)] = products_by_bacon2.get((rc1, rc2), RCGraphModule()) + coeff * rc
    #     if rc1.is_principal and rc2.is_principal:
    #         coproducts_by_bacon[rc.perm] = coproducts_by_bacon.get(rc.perm, TensorModule()) + coeff * TensorModule({RCGraphTensor(rc1,rc2): 1})

    for (rc, ((perm1, _), (perm2, _))), coeff in solution_module3.items():
        #perm1, perm2 = key1[0], key2[0]
        # perm = rc.perm
        if len(perm1) > n or len(perm2) > n or len(rc.perm) > n:
            continue
        products_rc[(perm1, perm2)] = products_rc.get((perm1, perm2), RCGraphModule()) + coeff * rc
        #products[(rc1.perm, rc2.perm)] = products_rc.get((rc1.perm, rc2.perm), 0) + coeff * Sx(rc.polyvalue(x))
        # assert coeff >= 0 NOPE
        # coprods_interim[rc.perm] = coprods_interim.get(rc.perm, TensorModule()) + coeff * TensorModule.ext_multiply(1*rc1, ASx(perm2, len(rc2)))
        # coprods_interim_rc[rc.perm] = coprods_interim.get(rc.perm, TensorModule()) + coeff * TensorModule.ext_multiply(Sx(rc1.perm), 1 * rc2)
        # coprods_length_vector[rc.length_vector()] = coprods_length_vector.get(rc.length_vector(), 0) + coeff * (FA@FA)((rc1.length_vector(), rc2.length_vector()))
        # coprods_length_vector2[rc.length_vector()] = coprods_length_vector2.get(rc.length_vector(), 0) + coeff * (ASx@ASx)(((rc1.perm,len(rc1)), (rc2.perm,len(rc2))))
        # if rc1.is_principal and rc2.is_principal:
        #     products_by_weight[(rc1.perm, rc2.perm)] = products_by_weight.get((rc1.perm, rc2.perm), RCGraphModule()) + coeff * rc

    # for (perm1, perm2), module in products_by_weight.items():
    #     print(perm1.trimcode, perm2.trimcode)
    #     print(module)

    
    # for (perm1, perm2), module in products.items():
    #     for rc, coeff in module.items():
    #         weight_sum[rc.length_vector()] = weight_sum.get(rc.length_vector(), 0) + coeff

    #process products_rc2
    

    

        
    # for tup, module in products_by_weight.items():
    #     for rc, coeff in module.items():
    #         if len(rc.perm) > n:
    #             continue
    #         assert rc.is_principal or coeff == 0, f"{coeff=} {rc.perm=} {rc.length_vector()=} {rc.perm.trimcode=}" # all but principal cancel
            #products_rc2[tup] = products_rc2.get(tup, 0) + coeff * Sx(rc.polyvalue(x))
        #products_rc[(rc1, rc2)] = products_rc.get((rc1, rc2), 0) + coeff * Sx(perm)

    assert all(v > 0 for v in weight_sum.values()), f"{weight_sum}" # total sums all >= 0

    # for length_vector, module in rc_sum.items():
    #     for rc, coeff in module.items():
    #         assert coeff == 0 or rc.is_principal, f"{coeff=} {rc.perm=} {rc.length_vector()=} {rc.perm.trimcode=}" # all but principal cancel
        #products_rc[(rc1, rc2)] = products_rc.get((rc1, rc2), 0) + coeff * Sx(perm)


    # for perm, module in coprods_interim.items():
    #     #assert all(v >=0 for v in module.values()), "NOPE"
    #     for (rc, key), coeff in module.items():
    #         coprods[perm] = coprods.get(perm, 0) + coeff * (ASx@ASx)(((rc.perm, len(rc)),key))


    # for perm, module in coprods_interim_rc.items():
    #     assert all(v >=0 for v in module.values()), "NOPE"
    #     for (perm1, rc2), coeff in module.items():
    #         coprods[perm] = coprods.get(perm, 0) + coeff * (ASx@ASx)(((perm1, n-1),(rc2.perm, len(rc2))))

    # for length_vector, sx_elem in coprods_length_vector2.items():
    #     print(sx_elem)
    #     #assert all(v >= 0 for v in fa_elem.values()), "NOPE"
    #     elem = FA(*length_vector).change_basis(SchubertBasis)
    #     for key, v in elem.items():
    #         assert v * coeff >= 0, "NOPE"
    #         coprods[key] = coprods.get(key, 0) + coeff * v * sx_elem

    # for length_vector, fa_elem in coprods_length_vector.items():
    #     print(fa_elem)
    #     #assert all(v >= 0 for v in fa_elem.values()), "NOPE"
    #     elem = FA(*length_vector).change_basis(SchubertBasis)
    #     for key, v in elem.items():
    #         assert v * coeff >= 0, "NOPE"
    #         coprods[key] = coprods.get(key, 0) + coeff * v * 1

    num_successes = 0

    # print("CHECK1!!!\n\n\n")

    # for perm, elem in coprods.items():
    #     print(f"{perm.trimcode}")
    #     print(elem)
    #     check = ASx(perm, n - 1).coproduct()
    #     print(check)
    #     diff = elem - check
    #     print(diff)
    #     assert all(v == 0 for v in diff.values()), f"Failed check on {perm}, diff = {diff}"
    #     num_successes += 1


    # print("CHECK2!!!\n\n\n")

    for (perm1, perm2), elem in products.items():
        if product_too_big(perm1, perm2, n):
            continue
        print(f"{perm1.trimcode}, {perm2.trimcode}")
        check = Sx(perm1)*Sx(perm2)
        sumup = 0
        
        try:
            assert elem == check
        except AssertionError as e:
            print(f"Fail")
            print(f"{elem=}")
            print(f"{check=}")
            continue
        print(f"Success {elem}")
        num_successes += 1

    #print("CHECK!!!\n\n\n")

    for (perm1, perm2), elem in products_rc.items():
        if product_too_big(perm1, perm2, n):
            continue
        print(f"{perm1.trimcode}, {perm2.trimcode}")
        check = Sx(perm1)*Sx(perm2)
        sumup = 0
        for rc, coeff in elem.items():
            try:
                assert check.get(rc.perm, 0) == coeff
            except AssertionError as e:
                print(f"Fail on {rc.perm.trimcode} {coeff=} {check.get(rc.perm, 0)=} {check=}")
                print(elem)
                continue
            if rc.is_principal:
                sumup += coeff * Sx(rc.perm)
        assert sumup == check
        print(f"Success {sumup}")
        num_successes += 1

    # print("CHECK4!!!\n\n\n")

    # for (perm1, perm2), elem in products_rc2.items():
    #     if product_too_big(perm1, perm2, n):
    #         continue
    #     print(f"{perm1.trimcode}, {perm2.trimcode}")
    #     check = Sx(perm1)*Sx(perm2)
    #     sumup = 0
    #     for rc, coeff in elem.items():
    #         # diff = elem - check
    #         # print(diff)
    #         assert check.get(rc.perm, 0) == coeff
    #         if rc.is_principal:
    #             sumup += coeff * Sx(rc.perm)
    #     assert sumup == check
    #     print(f"Success {sumup}")
    #     num_successes += 1

    # print("CHECK5!!!\n\n\n")

    # for (perm1, perm2), elem in products_by_weight.items():
    #     if product_too_big(perm1, perm2, n):
    #         continue
    #     print(f"{perm1.trimcode}, {perm2.trimcode}")
    #     check = Sx(perm1)*Sx(perm2)
    #     sumup = 0
    #     for rc, coeff in elem.items():
    #         # diff = elem - check
    #         # print(diff)
            
            
    #         if rc.is_principal:
    #             # print("princpy")
    #             # print(coeff)
    #             # assert check.get(rc.perm, 0) == coeff
    #             # #print(check[rc])
    #             # # md = FA(*rc.length_vector()).coproduct() * unit_tensor_rc_module
    #             # # addup = 0
    #             # # for (rc1, rc2), coeff0 in md.items():
    #             # #     if rc1.perm == perm1 and rc2.perm == perm2:
    #             # #         assert coeff0 == 1
    #             # #         print(f"adding {rc1.length_vector()=} {rc2.length_vector()=}")
    #             # #         addup += coeff0
    #             # # assert addup == coeff, f"{perm1=} {perm2=} {rc.perm=} {rc.length_vector()} {coeff=} {addup=}"
    #             sumup += coeff * Sx(rc.perm)

    #     assert sumup == check
    #     print(f"Success {sumup}")
    #     num_successes += 1

    
    # # for (rc1, rc2), elem in products_by_bacon.items():
    # #     # if not rc1.is_principal or not rc2.is_principal:
    # #     #     continue
    # #     perm1, perm2 = rc1.perm, rc2.perm
    # #     if product_too_big(perm1, perm2, n):
    # #         continue
    # #     print(f"{perm1.trimcode}, {perm2.trimcode}")
    # #     check = Sx(perm1)*Sx(perm2)
    #     #sumup = 0
    #     sumup = elem
    #     # for rc, coeff in elem.items():
    #     #     # diff = elem - check
    #     #     sumup += coeff * Sx(rc.perm)

    #     try:
    #         assert sumup == check, f"{sumup=} {check=}"
    #     except AssertionError as e:
    #         print(f"Fail these rcs")
    #         print(f"{RCGraphTensor(rc1,rc2)} {elem=} {check=}")
    #         continue
            
    #     print(f"Success {perm1.trimcode,perm2.trimcode} {sumup}")
    #     num_successes += 1

    # for (rc1, rc2), elem in products_by_bacon2.items():
    #     # if not rc1.is_principal or not rc2.is_principal:
    #     #     continue
    #     perm1, perm2 = rc1.perm, rc2.perm
    #     if product_too_big(perm1, perm2, n):
    #         continue
    #     print(f"{perm1.trimcode}, {perm2.trimcode}")
    #     check = Sx(perm1)*Sx(perm2)
    #     sumup = 0
    #     for rc, coeff in elem.items():
    #         # diff = elem - check
    #         sumup += coeff * Sx(rc.perm)
        
    #     try:
    #         assert sumup == check, f"{sumup=} {check=}"
    #     except AssertionError as e:
    #         print(f"fail {sumup=} {check=}")
    #         print(RCGraphTensor(rc1,rc2))
    #         continue
            
    #     print(f"Success {perm1.trimcode,perm2.trimcode} {sumup}")
    #     print(RCGraphTensor(rc1,rc2))
    #     num_successes += 1


    # for perm, elem in coproducts_by_bacon.items():
    #     # if not rc1.is_principal or not rc2.is_principal:
    #     #     continue
    #     #perm = rc.perm
    #     # perm1, perm2 = rc1.perm, rc2.perm
    #     # if product_too_big(perm1, perm2, n):
    #     #     continue
    #     print(f"{perm.trimcode}")
    #     check = ASx(perm, n - 1).coproduct()
    #     sumup = 0
    #     for (rc1, rc2), coeff in elem.items():
    #         # diff = elem - check
    #         if rc1.is_principal and rc2.is_principal:
    #             sumup += coeff * (ASx@ASx)(((rc1.perm, len(rc1)),(rc2.perm, len(rc2))))
        
    #     try:
    #         assert all(v==0 for v in (sumup - check).values()), f"{sumup=} {check=}"
    #     except AssertionError as e:
    #         print(f"fail")
    #         print(f"{sumup=}")
    #         print(f"{check=}")
    #         print(elem)
    #         continue
    #     print(f"Success {perm.trimcode}")
    #     print(elem)
    #     # print(f"Success {perm1.trimcode,perm2.trimcode} {sumup}")
    #     # print(RCGraphTensor(rc1,rc2))
    #     num_successes += 1


    # for (rc1, rc2), rc3 in products_by_bacon.items():
    #     if product_too_big(perm1, perm2, n):
    #         continue
    #     print(f"{perm1.trimcode}, {perm2.trimcode}")
    #     check = Sx(perm1)*Sx(perm2)
    #     sumup = 0
    #     for rc, coeff in elem.items():
    #         # diff = elem - check
    #         # print(diff)
            
            
    #         if rc.is_principal:
    #             # print("princpy")
    #             # print(coeff)
    #             # assert check.get(rc.perm, 0) == coeff
    #             # #print(check[rc])
    #             # # md = FA(*rc.length_vector()).coproduct() * unit_tensor_rc_module
    #             # # addup = 0
    #             # # for (rc1, rc2), coeff0 in md.items():
    #             # #     if rc1.perm == perm1 and rc2.perm == perm2:
    #             # #         assert coeff0 == 1
    #             # #         print(f"adding {rc1.length_vector()=} {rc2.length_vector()=}")
    #             # #         addup += coeff0
    #             # # assert addup == coeff, f"{perm1=} {perm2=} {rc.perm=} {rc.length_vector()} {coeff=} {addup=}"
    #             sumup += coeff * Sx(rc.perm)

    #     assert sumup == check
    #     print(f"Success {sumup}")
    #     num_successes += 1
    
    # print("CHECK4!!!\n\n\n")
    # for (perm1, perm2), elem in sumup.items():
    #     print(f"{perm1.trimcode}, {perm2.trimcode}")
    #     print(elem)
    #     check = Sx(perm1)*Sx(perm2)
    #     assert elem == checks
    #     print(f"Success {perm1.trimcode,perm2.trimcode}")
    #     num_successes += 1
    # for (perm1, perm2), elem in products_rc2.items():
    #     if product_too_big(perm1, perm2, n):
    #         continue
    #     print(f"{perm1.trimcode}, {perm2.trimcode}")
    #     check = Sx(perm1)*Sx(perm2)
    #     sumup = elem
    #     # for rc, coeff in elem.items():
    #     #     # diff = elem - check
    #     #     # print(diff)
    #     #     if rc.is_principal:
    #     #         assert check.get(rc.perm, 0) == coeff
    #     #     else:
    #     #         assert coeff == 0
    #     #     if rc.is_principal:
    #     #         sumup += coeff * Sx(rc.perm)
    #     assert sumup == check, f"{sumup=} {check=}"
    #     print(f"Success {sumup}")
    #     num_successes += 1

    # for (rc1, rc2), elem in products_rc.items():
    #     perm1, perm2 = rc1.perm, rc2.perm
    #     if product_too_big(perm1, perm2, n):
    #         continue
    #     print(f"{perm1.trimcode}, {perm2.trimcode}")
    #     print(elem)
    #     check = Sx(perm1)*Sx(perm2)
    #     print(check)
    #     # diff = elem - check
    #     # print(diff)
    #     assert elem == check
    #     num_successes += 1

    

    #assert num_successes == math.factorial(n), f"Only {num_successes} successes, not the full group"

if __name__ == "__main__":
    main()
