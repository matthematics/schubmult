import sys

from schubmult import ASx, Permutation, uncode
from schubmult.abc import x, y, z
from schubmult.rings import FA, FreeAlgebra, FreeAlgebraBasis, MonomialBasis, NilHeckeRing, PolynomialAlgebra, SchubertBasis, SingleSchubertRing, Sx, TensorRing
from schubmult.rings.rc_graph_module import RCGraph, RCGraphModule, TensorModule
from schubmult.symbolic import S, expand, expand_seq
from schubmult.utils.perm_utils import artin_sequences


def vector_sum(v1, v2):
    return tuple([a + b for a, b in zip(v1, v2)])

def main():
    n = int(sys.argv[1])

    seqs = artin_sequences(n - 1)
    xring = PolynomialAlgebra(MonomialBasis(x, n - 1))

    yring = PolynomialAlgebra(MonomialBasis(y,n-1))
    zring = PolynomialAlgebra(MonomialBasis(z,n-1))

    unit_rc_module = RCGraphModule({RCGraph(): 1})
    result = TensorModule()

    # 100% positive!
    #ASx([]).ring @ Sx([]).ring
    test_addup = TensorModule()
    if n == 3:
        assert (0,1) in seqs
    for seq1 in seqs:
        for seq2 in seqs:
            seq = vector_sum(seq1, seq2)
            if seq not in seqs:
                continue

        # left_result = TensorModule()
    #     addup = TensorModule()
            rc_graph_module_term = FA(*seq) * unit_rc_module
        #     #rc_graph_module_term = RCGraphModule({rc: v for rc, v in rc_graph_module_term.items() if len(rc.perm) <= n})
            tmodule = TensorModule.ext_multiply(rc_graph_module_term, rc_graph_module_term)
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
            module2 = FA(*seq1) * unit_rc_module
            module3 = FA(*seq2) * unit_rc_module
        #             schub_elem1 = FA(*seq1).change_basis(SchubertBasis)
        #             schub_elem2 = FA(*seq2).change_basis(SchubertBasis)
        #             #module2 = RCGraphModule({rc: v for rc, v in module2.items() if len(rc.perm) <= n})
        #             #module3 = RCGraphModule({rc: v for rc, v in module3.items() if len(rc.perm) <= n})
            tmodule2 = TensorModule.ext_multiply(module2, module2)
            tmodule3 = TensorModule.ext_multiply(module3, module3)

            for (rc1_left, rc1_right), coeff1 in tmodule.items():
                for (rc2_left, rc2_right), coeff2 in tmodule2.items():
                    for (rc3_left, rc3_right), coeff3 in tmodule3.items():
                        if rc1_left == rc1_right and rc2_left == rc2_right and rc3_left == rc3_right:
                            if rc1_left.perm == uncode(rc1_left.length_vector()) and rc2_left.perm == uncode(rc2_left.length_vector()) and rc3_left.perm == uncode(rc3_left.length_vector()):
                                new_coeff = coeff2 * coeff3 * coeff1
                                result += new_coeff * TensorModule.ext_multiply(ASx(rc1_left.perm, len(rc1_left)),TensorModule.ext_multiply(ASx(rc2_left.perm, len(rc2_left)), ASx(rc3_left.perm, len(rc3_left))))

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
    failed = False
    perms = set()
    perm_pairs = set()
    perm_pairs_total = set()
    rc_addups_left = {}
    rc_addups_right = {}
    
    sum1 = {}
    sum2 = {}
    full_addups = {}
    # saw = False
    addups = {}
    for ((perm, _), ((left_coprod_perm, _), (right_coprod_perm, _))), value in result.items():
        # if any(len(permperm) > n for permperm, val in (Sx(left_coprod_perm) * Sx(right_coprod_perm)).items() if val != S.Zero):
        #     continue
        if len(perm) > n or len(left_coprod_perm) > n or len(right_coprod_perm) > n:
            continue
        perm_pairs_total.add((left_coprod_perm, right_coprod_perm))
        #addups[(left_coprod_perm, right_coprod_perm)] = addups.get((left_coprod_perm, right_coprod_perm), RCGraphModule()) + RCGraphModule({result_rc: value})
        # if perm != uncode(vector_sum(rc1.length_vector(), rc2.length_vector())):# or rc1.perm != uncode(rc1.length_vector()) or rc2.perm != uncode(rc2.length_vector()):
        #     continue
        product = Sx(left_coprod_perm) * Sx(right_coprod_perm)
        #RCGraphModule(dict.fromkeys(RCGraph.all_rc_graphs(result_rc.perm, n-1), value))
        
        testval = value - product.get(perm, 0)
        if testval != S.Zero:
            print(f"Failures for {(left_coprod_perm.trimcode, right_coprod_perm.trimcode)} {perm.trimcode} {value=} {product.get(perm, 0)=} {testval=}")
            print(perm)
            failed = True
        else:
            print(f"Success for {(left_coprod_perm.trimcode, right_coprod_perm.trimcode)}: Sx({left_coprod_perm.trimcode})*Sx({right_coprod_perm.trimcode})")
            print(perm)
            perms.add(perm)
            perm_pairs.add((left_coprod_perm, right_coprod_perm))
    # assert saw
    # side_addups = {}
    # side_addups2 = {}
    # for key in rc_addups_left:
    #     assert expand(rc_addups_left[key].polyvalue(x) -  rc_addups_right[key].polyvalue(x)) == S.Zero, f"Failure for {key}: {rc_addups_left[key].polyvalue(x)} != {rc_addups_right[key]}"
    #     print(f"{key[0].perm.trimcode,key[1].perm.trimcode,key[2].perm.trimcode}")
    #     print(rc_addups_left[key])
    #     print(rc_addups_right[key])
    #     print(full_addups[key])
    #     side_addups[(key[0].perm,key[1].perm)] = side_addups.get((key[0].perm,key[1].perm), RCGraphModule()) + rc_addups_left[key]
    #     side_addups2[key[2].perm] = side_addups2.get(key[2].perm, TensorModule()) + rc_addups_right[key]
    

    # keys = {(key[0].perm,key[1].perm,key[2].perm) for key in rc_addups_left}

    # for key in keys:
    #     print(f"FINAL FOR {key[0].trimcode,key[1].trimcode}")
    #     print(side_addups[(key[0],key[1])])
    #     print(f"FINAL BUBBLE FOR {key[2].trimcode}")
    #     print(side_addups2[key[2]])
    
    # for (perm1, perm2), module in addups.items():
    #     print(f"THIS IS {(perm1.trimcode, perm2.trimcode)}")
    #     print(module)
 

        

    # for (rc1,rc2), module in rc_addups.items():
    #     print(f"THIS IS {(rc1.trimcode, rc2.trimcode)}")
    #     print(module)

        

    # for perm in perm_pairs_total:
    #     if perm not in perm_pairs:
    #         assert len(uncode([a+b for a,b in zip(perm[0].trimcode, perm[1].trimcode)])) > n, f"Failure for {perm}"
            #print(result_dict[coprod_key])
    #print(new_result)
    # exit()
    # result_dict = {}
    # for key, value in result.items():
    #     perm1, perm2 = key[0]
    #     coprod_perm_pair = key[1]
    #     if perm1 == perm2:
    #         result_dict[coprod_perm_pair] = result_dict.get(coprod_perm_pair, 0) + value * Sx(perm1)

    # Permutation.print_as_code = True
    # result_dict = {}
    # failed = False
    # for key, value in result.items():
    #     #right_graph = key[0][1][1]
    #     right_nil_perm = key[0][1][1]
    #     #$right_nil_perm = key[0][1][1].perm
    #     #left_graph = key[0][0][0]
    #     left_nil_perm = key[0][0][0]
    #     right_schub_perm = key[0][1][0]
    #     coprod_perm_pair = key[1]
    #     left_schub_perm = key[0][0][1]
    #     if right_nil_perm == right_schub_perm and left_nil_perm == right_nil_perm:
            
    
    # for coprod_key in result_dict:
    #     ((left_coprod_perm, _), (right_coprod_perm, _)) = coprod_key
    #     if any(len(permperm) > n for permperm, val in (Sx(left_coprod_perm) * Sx(right_coprod_perm)).items() if val != S.Zero):
    #         continue
    #     testval = Sx(result_dict[coprod_key].polyvalue(x)) - Sx(left_coprod_perm) * Sx(right_coprod_perm)
    #     if testval.expand() != S.Zero:
    #         print(f"Failures for {(left_coprod_perm.trimcode, right_coprod_perm.trimcode)}: {testval=} {result_dict[coprod_key]=}")
    #         failed = True
    #     else:
    #         print(f"Success for {(left_coprod_perm.trimcode, right_coprod_perm.trimcode)}: Sx({left_coprod_perm.trimcode})*Sx({right_coprod_perm.trimcode})={result_dict[coprod_key]}")

    if not failed:
        print(f"YEAH!!! {len(perms)=} {len(perm_pairs)=} {len(perm_pairs_total)=}", file=sys.stderr)
        print(f"YEAH!!! {len(perms)=} {len(perm_pairs)=} {len(perm_pairs_total)=}", file=sys.stderr)


if __name__ == "__main__":
    main()
