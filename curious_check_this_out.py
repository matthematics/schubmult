import sys

from schubmult import ASx, Permutation, uncode
from schubmult.abc import x, y, z
from schubmult.rings import FA, FreeAlgebra, FreeAlgebraBasis, MonomialBasis, NilHeckeRing, PolynomialAlgebra, SchubertBasis, SingleSchubertRing, Sx, TensorRing, WordBasis
from schubmult.rings.rc_graph_module import RCGraph, RCGraphModule, RCGraphTensor, TensorModule, all_fa_degree
from schubmult.symbolic import S, expand, expand_seq
from schubmult.utils.perm_utils import artin_sequences


def vector_sum(v1, v2):
    return tuple([a + b for a, b in zip(v1, v2)])

def rc_ify(rc_graph):
            return TensorModule.ext_multiply(ASx(rc_graph.perm, len(rc_graph.length_vector())), RCGraphModule({rc_graph: 1}))
        
def elem_ify(rc_graph):
    return TensorModule.ext_multiply(ASx(rc_graph.perm, len(rc_graph)), Sx(RCGraphModule({rc_graph: 1}).polyvalue(x)))

def pos_mod(seq, n):
    unit_module = RCGraphModule({RCGraph(): 1})
    ret = FA(*seq) * unit_module
    return RCGraphModule({rc: v for rc, v in ret.items() if len(rc.perm) <= n})

def seq_schubert(seq, n):
    unit_module = RCGraphModule({RCGraph(): 1})
    seq_module = FA(*seq) * unit_module
    # for seq0 in artin_sequences(len(seq)):
    seq0 = seq
    seq0_module = FA(*seq0) * unit_module
    ret_module = TensorModule()
    basis_change = {}
    for rc in seq0_module:
        base = ASx(rc.perm, len(rc)).change_basis(WordBasis)
        basis_change[rc.perm] = base.get(seq, 0)
    for rc, coeff in seq0_module.items():
        if len(rc.perm) > n:
            continue
        for rc2, coeff2 in seq_module.items():
            if len(rc2.perm) > n:
                continue
            if basis_change[rc.perm] != 0:
                ret_module += TensorModule({RCGraphTensor(rc.perm, rc2): coeff * coeff2 *  basis_change[rc.perm]})
    return ret_module


# def seq_schabert(seq, n):
#     unit_module = RCGraphModule({RCGraph(): 1})
#     seq_module = FA(*seq) * unit_module
#     ret_module = TensorModule()
#     basis_change = {}
#     # for rc in seq_module:
#     #     base = ASx(rc.perm, len(rc)).change_basis(WordBasis)
#     #     basis_change[rc.perm] = base.get(rc.length_vector(), 0)
#     for rc, coeff in seq_module.items():
#         if len(rc.perm) > n:
#             continue
#         for rc2, coeff2 in seq_module.items():
#             if len(rc2.perm) > n:
#                 continue
#             #if basis_change[rc.perm] != 0:
#             if rc2 == rc:
#                 ret_module += TensorModule({RCGraphTensor(rc, rc2): coeff * coeff2})
#     return ret_module

def main():
    n = int(sys.argv[1])

    seqs = artin_sequences(n - 1)
    deg = (n * (n - 1)) // 2
    

    unit_rc_module = RCGraphModule({RCGraph(): 1})
    result = TensorModule()

    # 100% positive!
    #ASx([]).ring @ Sx([]).ring
    test_addup = TensorModule()
    
    test_module = 0
    #seqs = all_fa_degree(degree, n - 1)
    for seq in seqs:
        for rc, coeff in pos_mod(seq, n).items():
            test_module += coeff * (ASx([]).ring@Sx([]).ring).from_dict(elem_ify(rc))
        smod = seq_schubert(seq, n)
        #print(f"{seq=}")
        #print(smod)
        #test_module += seq_schubert(seq, n)
        #print(test_module)
        # if smod == 0:
        #     continue

    #exit()
        # for seq2 in seqs:
        #     seq = vector_sum(seq1, seq2)
        #     if seq not in seqs:
        #         continue
        coprod_base = FA(*seq).coproduct()
        #coprod = FreeAlgebraBasis.change_tensor_basis(coprod_base, SchubertBasis, SchubertBasis)
    # left_result = TensorModule()
#     addup = TensorModule()
            #rc_graph_module_term = FA(*seq) * unit_rc_module
    #     #rc_graph_module_term = RCGraphModule({rc: v for rc, v in rc_graph_module_term.items() if len(rc.perm) <= n})
        #tmodule = TensorModule.ext_multiply(rc_graph_module_term, rc_graph_module_term)
    #         #if rc1 == rc2:
    # #                addup += TensorModule.ext_multiply(ASx(rc1.perm,len(rc1)), RCGraphModule({rc2: coeff}))

    #     coprod = FA(*seq).coproduct()
    #     #test_addup += addup
    #     #left_result = TensorModule.ext_multiply(addup, FreeAlgebraBasis.change_tensor_basis(coprod, SchubertBasis, SchubertBasis))
    #     left_result = addup
    #     # # zr_elem = zring(seq)
    #     # # xr_elem = xring(seq)
    #     # # module1 = TensorModule.ext_multiply(rc_graph_module_term, zr_elem)
    #     # # module2 = TensorModule.ext_multiply(xr_elem, rc_graph_module_term)
    #     # module1 = TensorModule.ext_multiply(rc_graph_module_term, rc_graph_module_term)
    #     # # module2 = TensorModule.ext_multiply(rc_graph_module_term, rc_graph_module_term)
    #     right_result = TensorModule()
    #     for (rc1, rc2), coeff00 in tmodule.items():
    #         for (seq1, seq2), val in coprod.items():

        # for (rc1_left, rc1_right), coeff1 in tmodule.items():
        #     if len(rc1_left.perm) > n or len(rc1_right.perm) > n:
        #         continue
        #     if rc1_left != rc1_right:
        #         continue    
        for rc, coeff in pos_mod(seq, n).items():
            for (seq1, seq2), val in coprod_base.items():
            #         #result += val * TensorModule.ext_multiply(smod, TensorModule.ext_multiply(seq_schubert(seq1, n), seq_schubert(seq2, n)))
            #         #result += val * TensorModule.ext_multiply(TensorModule.ext_multiply(pos_mod(seq, n), smod), TensorModule.ext_multiply(TensorModule.ext_multiply(pos_mod(seq1,n),seq_schubert(seq1, n)), TensorModule.ext_multiply(pos_mod(seq2,n),seq_schubert(seq2, n))))
            #         #result += val * TensorModule.ext_multiply(TensorModule.ext_multiply(pos_mod(seq,n), smod), TensorModule.ext_multiply(TensorModule.ext_multiply(pos_mod(seq1,n),seq_schubert(seq1, n)), TensorModule.ext_multiply(pos_mod(seq2,n),seq_schubert(seq2, n))))
            #         #result += val * TensorModule.ext_multiply(smod, TensorModule.ext_multiply(seq_schubert(seq1,n),seq_schubert(seq2,n)))
            #         result += coeff * val * TensorModule.ext_multiply(seq_schubert(rc), TensorModule.ext_multiply(seq_schubert(seq1,n),seq_schubert(seq2,n)))
                for rc1, coeff1 in pos_mod(seq1, n).items():
                    for rc2, coeff2 in pos_mod(seq2, n).items():
                        # if len(rc1.perm) > n or len(rc2.perm) > n:
                        #     continue
                        # if rc1 != rc2:
                        #     continue
                        #test_addup += val * coeff * coeff2 * TensorModule.ext_multiply(smod, TensorModule.ext_multiply(elem_ify(rc1),elem_ify(rc2)))
                        #test_module += val * coeff * coeff2 * TensorModule.ext_multiply(smod, TensorModule.ext_multiply(rc_ify(rc1),rc_ify(rc2)))
                        #result += coeff * coeff2 * val * TensorModule.ext_multiply(smod, TensorModule.ext_multiply(elem_ify(rc1),elem_ify(rc2)))
                        result += coeff * coeff1 * coeff2 * val * TensorModule.ext_multiply(elem_ify(rc), TensorModule.ext_multiply(elem_ify(rc1),elem_ify(rc2)))
    # print("result=")
    # print(result)
    print(test_module)
    failed = False
    perm_pairs_total = set()
    perm_trips = set()
    perms = set()
    addup = {}
    #for ((sumrc1,(rc1_left, rc1_right)), ((sumrc2, (rc2_left, rc2_right)),(sumrc3, (rc3_left,rc3_right)))), coeff in result.items():
    #for ((sumrc1,(rc1_left, rc1_right)), (rc2_left, rc3_left)), coeff in result.items():
    for ((a, rc1_left), ((b, rc2_left), (c, rc3_left))), coeff in result.items():
        if a[0] == rc1_left and b[0] == rc2_left and c[0] == rc3_left:
            print(rc1_left)
            product = Sx(rc2_left) * Sx(rc3_left)
            assert product.get(rc1_left, 0) == coeff, f"Failure for left {(rc1_left.trimcode, rc2_left.trimcode, rc3_left.trimcode)}: {product.get(rc1_left, 0)=} {coeff=} {rc1_left.trimcode=}"
            print(f"Success {rc1_left.trimcode, rc2_left.trimcode, rc3_left.trimcode}: Sx({rc2_left.trimcode})*Sx({rc3_left.trimcode})=...+{coeff * Sx(rc1_left.trimcode)}")
            perm_pairs_total.add((rc2_left, rc3_left))
            perms.add(rc1_left)
            perm_trips.add((rc1_left, rc2_left, rc3_left))
            st = [rc2_left, rc3_left]
        # if sumrc2.perm == rc2_left.perm and sumrc3.perm == rc3_left.perm:
            addup[tuple(st)] = addup.get(tuple(st), 0) + coeff * Sx(rc1_left)

    for (perm1, perm2), val in addup.items():
        product = Sx(perm1) * Sx(perm2)
        try:
            assert product == val
            print(f"Success for {(perm1.trimcode, perm2.trimcode)}: {val=}")
        except AssertionError:
            assert any(len(perm) > n for perm in product), f"Failure for {(perm1.trimcode, perm2.trimcode)}: {product=} {val=}"
    
    # for (perm1, perm2, perm3) in perm_trips:
    #     assert all([c + b == a for a, b, c in zip(perm1.trimcode, perm2.trimcode, perm3.trimcode)]), f"Failure for {(perm1.trimcode, perm2.trimcode, perm3.trimcode)}"
    # test_dict = {}
    # for ((elem, rc), ((f_elem1, elem1), (f_elem2, elem2))), coeff in result.items():
    #     assert elem[0] == rc.perm
    #     assert elem1 == f_elem1[0]
    #     assert elem2 == f_elem2[0], f"{elem2=}, {f_elem2=} {coeff=}"
    #     #if uncode(rc.length_vector()) == elem[0]:
    #     test_dict[(elem1, elem2)] = test_dict.get((elem1, elem2), 0) + coeff * Sx(rc.polyvalue(x))

    # for (perm1, perm2), module in test_dict.items():
    #     assert Sx(perm1) * Sx(perm2) == module, f"Failure for {(perm1.trimcode, perm2.trimcode)}: {Sx(perm1) * Sx(perm2)=} {module=}"
    # print("Passed the first test!")


    # old_result = result
    # result = TensorModule()
    # for ((elem, rc), ((elem1, rc1), (elem2, rc2))), coeff in old_result.items():
    #     if len(rc.perm) > n or len(rc1.perm) > n or len(rc2.perm) > n:
    #         continue
    #     if rc.perm == uncode(rc.length_vector()):
    #         result += TensorModule.ext_multiply(coeff * Sx(rc.perm),TensorModule.ext_multiply(rc_ify(rc1),rc_ify(rc2)))

    # old_result = result
    # result = TensorModule()
    # for (elem, ((elem1, rc1), (elem2, rc2))), coeff in old_result.items():
    #     if rc1.perm == uncode(rc1.length_vector()):
    #         result += TensorModule.ext_multiply(coeff * Sx(elem),TensorModule.ext_multiply(Sx(rc1.perm),rc_ify(rc2)))

    # old_result = result
    # result = TensorModule()
    # for (elem, (elem1, (elem2, rc2))), coeff in old_result.items():
    #     if rc2.perm == uncode(rc2.length_vector()):
    #         result += TensorModule.ext_multiply(coeff * Sx(elem),TensorModule.ext_multiply(Sx(elem1),Sx(rc2.perm)))
    #             full_tmodule = TensorModule.ext_multiply(tmodule2,tmodule3)
    #             addup = TensorModule()
    #             for ((rc1_left, rc2_left),(rc1_right,rc2_right)), coeff0 in full_tmodule.items():
    #                 coeff = coeff0 * coeff00
    #                 if rc1 == rc2 and rc1_left == rc2_left and rc1_right == rc2_right:
    #                     #print(f'{rc1_left.perm.trimcode, rc1_right.perm.trimcode}, coeff={coeff}, val={val}, coeff1={coeff1}, coeff2={coeff2}, seq1={seq1}, seq2={seq2}')
    #                     # print(f"{schub_elem1=}, {schub_elem2=}")
    #                     new_coeff = coeff * schub_elem1.get((rc1_left.perm, len(rc1_left)), 0) * schub_elem2.get((rc1_right.perm, len(rc1_right)), 0)
    #                     #print(f"{new_coeff=}")
    #                     #addup += TensorModule.ext_multiply(TensorModule.ext_multiply(RCGraphModule({rc1_left: new_coeff}), RCGraphModule({rc1_right: 1})), 
    #                     #                                    TensorModule.ext_multiply(ASx(rc1_left.perm,len(rc1_left)), ASx(rc1_right.perm,len(rc1_right))))
    #                     result += TensorModule.ext_multiply(TensorModule.ext_multiply(ASx(rc1.perm,len(rc1)), RCGraphModule({rc2: new_coeff})),
    #                                                         TensorModule.ext_multiply(TensorModule.ext_multiply(RCGraphModule({rc1_left: coeff}), RCGraphModule({rc1_right: 1})),
    #                                                         TensorModule.ext_multiply(ASx(rc1_left.perm,len(rc1_left)), ASx(rc1_right.perm,len(rc1_right)))))
                #right_result += addup
        #                                                                                             TensorModule.ext_multiply(module3,module3)),
        #                            FreeAlgebraBasis.change_tensor_basis((FA @ FA)((seq1,seq2)),SchubertBasis,SchubertBasis)))
        #result += TensorModule.ext_multiply(left_result,right_result)
    # buildup_mod = {}
    # for ((perm, _), rc), coeff in test_addup.items():
    #     assert Sx(rc.polyvalue(x)) == Sx(perm), f"Failure for {perm.trimcode}"
    #     assert all(rc.perm == perm for rc in test_addup[rc].keys()), f"Failure for {perm.trimcode}"
    # left_graph == right_graph
    #result2 = result

    # result_dict = {}
    # for ((perm, _),rc), coeff in addup.items():
    #     result_dict[perm] = result_dict.get(perm, RCGraphModule()) + RCGraphModule({rc: coeff})

    # for perm in result_dict:
    #     assert Sx(result_dict[perm].polyvalue(x)) == Sx(perm), f"Failure for {perm.trimcode}"
    #     assert all(rc.perm == perm for rc in result_dict[perm].keys()), f"Failure for {perm.trimcode}"
    # exit()
    # new_result = TensorModule()
    
    # perms = set()
    # perm_pairs = set()
    # perm_pairs_total = set()
    # rc_addups_left = {}
    # rc_addups_right = {}
    
    # sum1 = {}
    # sum2 = {}
    # full_addups = {}
    # # saw = False
    # addups = {}
    # #for (perm, (rc1, rc2)), value in result.items():
    # for (result_rc, (((left_coprod_perm, _), perm1), ((right_coprod_perm, _), perm2))), value in result.items():
    #     perm = result_rc.perm
    #     #if any(len(permperm) > n for permperm, val in (Sx(rc1.perm) * Sx(rc2.perm)).items() if val != S.Zero):
    #     if any(len(permperm) > n for permperm, val in (Sx(left_coprod_perm) * Sx(right_coprod_perm)).items() if val != S.Zero):
    #         continue
    #     # left_coprod_perm, right_coprod_perm = rc1.perm, rc2.perm
    #     # if len(perm) > n or len(left_coprod_perm) > n or len(right_coprod_perm) > n:
    #     #     continue
        
    #     perm_pairs_total.add((left_coprod_perm, right_coprod_perm))
    #     #addups[(left_coprod_perm, right_coprod_perm)] = addups.get((left_coprod_perm, right_coprod_perm), RCGraphModule()) + RCGraphModule({result_rc: value})
    #     if result_rc.perm != uncode(result_rc.length_vector()):
    #         continue
    #     product = Sx(left_coprod_perm) * Sx(right_coprod_perm)
    #     #RCGraphModule(dict.fromkeys(RCGraph.all_rc_graphs(result_rc.perm, n-1), value))
    #     addups[(left_coprod_perm, right_coprod_perm)] = addups.get((left_coprod_perm, right_coprod_perm), 0) + value * Sx(perm)
    #     testval = value - product.get(perm, 0)
    #     if testval != S.Zero:
    #         print(f"Failures for {(left_coprod_perm.trimcode, right_coprod_perm.trimcode)} {perm.trimcode} {value=} {product.get(perm, 0)=} {testval=}")
    #         print(perm)
    #         failed = True
    #     else:
    #         print(f"Success for {(left_coprod_perm.trimcode, right_coprod_perm.trimcode)}: Sx({left_coprod_perm.trimcode})*Sx({right_coprod_perm.trimcode})")
    #         print(perm)
    #         perms.add(perm)
    #         perm_pairs.add((left_coprod_perm, right_coprod_perm))

    # for (perm1, perm2), module in addups.items():
    #     assert Sx(perm1) * Sx(perm2) == module, f"Failure for {(perm1.trimcode, perm2.trimcode)}: {Sx(perm1) * Sx(perm2)=} {module=}"
    # # assert saw
    # # side_addups = {}
    # # side_addups2 = {}
    # # for key in rc_addups_left:
    # #     assert expand(rc_addups_left[key].polyvalue(x) -  rc_addups_right[key].polyvalue(x)) == S.Zero, f"Failure for {key}: {rc_addups_left[key].polyvalue(x)} != {rc_addups_right[key]}"
    # #     print(f"{key[0].perm.trimcode,key[1].perm.trimcode,key[2].perm.trimcode}")
    # #     print(rc_addups_left[key])
    # #     print(rc_addups_right[key])
    # #     print(full_addups[key])
    # #     side_addups[(key[0].perm,key[1].perm)] = side_addups.get((key[0].perm,key[1].perm), RCGraphModule()) + rc_addups_left[key]
    # #     side_addups2[key[2].perm] = side_addups2.get(key[2].perm, TensorModule()) + rc_addups_right[key]
    

    # # keys = {(key[0].perm,key[1].perm,key[2].perm) for key in rc_addups_left}

    # # for key in keys:
    # #     print(f"FINAL FOR {key[0].trimcode,key[1].trimcode}")
    # #     print(side_addups[(key[0],key[1])])
    # #     print(f"FINAL BUBBLE FOR {key[2].trimcode}")
    # #     print(side_addups2[key[2]])
    
    # # for (perm1, perm2), module in addups.items():
    # #     print(f"THIS IS {(perm1.trimcode, perm2.trimcode)}")
    # #     print(module)
 

        

    # # for (rc1,rc2), module in rc_addups.items():
    # #     print(f"THIS IS {(rc1.trimcode, rc2.trimcode)}")
    # #     print(module)

        

    # # for perm in perm_pairs_total:
    # #     if perm not in perm_pairs:
    # #         assert len(uncode([a+b for a,b in zip(perm[0].trimcode, perm[1].trimcode)])) > n, f"Failure for {perm}"
    #         #print(result_dict[coprod_key])
    # #print(new_result)
    # # exit()
    # # result_dict = {}
    # # for key, value in result.items():
    # #     perm1, perm2 = key[0]
    # #     coprod_perm_pair = key[1]
    # #     if perm1 == perm2:
    # #         result_dict[coprod_perm_pair] = result_dict.get(coprod_perm_pair, 0) + value * Sx(perm1)

    # # Permutation.print_as_code = True
    # # result_dict = {}
    # # failed = False
    # # for key, value in result.items():
    # #     #right_graph = key[0][1][1]
    # #     right_nil_perm = key[0][1][1]
    # #     #$right_nil_perm = key[0][1][1].perm
    # #     #left_graph = key[0][0][0]
    # #     left_nil_perm = key[0][0][0]
    # #     right_schub_perm = key[0][1][0]
    # #     coprod_perm_pair = key[1]
    # #     left_schub_perm = key[0][0][1]
    # #     if right_nil_perm == right_schub_perm and left_nil_perm == right_nil_perm:
            
    
    # # for coprod_key in result_dict:
    # #     ((left_coprod_perm, _), (right_coprod_perm, _)) = coprod_key
    # #     if any(len(permperm) > n for permperm, val in (Sx(left_coprod_perm) * Sx(right_coprod_perm)).items() if val != S.Zero):
    # #         continue
    # #     testval = Sx(result_dict[coprod_key].polyvalue(x)) - Sx(left_coprod_perm) * Sx(right_coprod_perm)
    # #     if testval.expand() != S.Zero:
    # #         print(f"Failures for {(left_coprod_perm.trimcode, right_coprod_perm.trimcode)}: {testval=} {result_dict[coprod_key]=}")
    # #         failed = True
    # #     else:
    # #         print(f"Success for {(left_coprod_perm.trimcode, right_coprod_perm.trimcode)}: Sx({left_coprod_perm.trimcode})*Sx({right_coprod_perm.trimcode})={result_dict[coprod_key]}")

    if not failed:
        print(f"YEAH!!! {len(perms)=} {len(perm_pairs_total)=} {len(perm_trips)=}", file=sys.stderr)
        print(f"YEAH!!! {len(perms)=} {len(perm_pairs_total)=} {len(perm_trips)=}", file=sys.stderr)


if __name__ == "__main__":
    main()

