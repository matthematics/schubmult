import sys

from schubmult import ASx, Permutation, uncode
from schubmult.abc import x
from schubmult.rings import FA, Sx, WordBasis
from schubmult.rings.free_algebra_basis import FreeAlgebraBasis, SchubertBasis
from schubmult.rings.rc_graph_module import RCGraph, RCGraphModule, RCGraphTensor, TensorModule, all_fa_degree
from schubmult.symbolic import expand_seq
from schubmult.utils.perm_utils import artin_sequences


def vector_sum(v1, v2):
    return tuple([a + b for a, b in zip(v1, v2)])


def product_too_big(perm1, perm2, n):
    return any(len(perm)>n for perm in (Sx(perm1) * Sx(perm2)).keys())

def main():
    n = int(sys.argv[1])

    unit_rc_module = RCGraphModule({RCGraph(): 1})

    # 100% positive!

    unit_tensor_rc_module = TensorModule.ext_multiply(unit_rc_module, unit_rc_module)

    solution_module = TensorModule()
    solution_module2 = TensorModule()
    solution_module3 = TensorModule()


    # this is an inner product of RCs
    perms = Permutation.all_permutations(n)
    # Act ASx, principals cancel
    aseqs = artin_sequences(n-1)
    degree = (n*(n-1))//2
    seqs = set()
    for deg in range(degree+1):
        seqs.update(all_fa_degree(deg, n-1))
    # for seq in aseqs:
    # #for perm in aseqs:

    #     pish_mod = FA(*seq) * unit_rc_module
    #     #pish_mod = FA(*seq) * unit_rc_module
    #     for rc, bob in pish_mod.items():
    #         if len(rc.perm) > n:
    #             continue
    #         #solution_module += coeff*TensorModule.ext_multiply(pish_mod,ASx(rc.perm, n-1).coproduct())
    #         # fist = ASx(rc.perm, n-1).coproduct()*unit_tensor_rc_module
    #         # for (rc1, rc2), coeff2 in fist.items():
    #         #     if vector_sum(rc1.length_vector(), rc2.length_vector()) == seq:
    #         #         #print(f"Adding {(rc, (rc1, rc2))} with coeff {bob*coeff*coeff2}")
    #         #solution_module2 += bob * TensorModule.ext_multiply(pish_mod,ASx(rc.perm,n-1).coproduct()*unit_tensor_rc_module)
    #         # did the whole pish mod!!!!
    #         solution_module2 += TensorModule.ext_multiply(bob * rc,ASx(rc.perm,len(rc)).coproduct()*unit_tensor_rc_module)
    solution_module = TensorModule()
    pickle_module = TensorModule()
    # for perm in perms:
    #     acter = ASx(perm, n-1).change_basis(WordBasis)
    #     graphs = RCGraph.all_rc_graphs(perm, n - 1)
    #     graph_module = RCGraphModule(dict.fromkeys(graphs, 1))
    #     for word, coeff in acter.items():
    #         solution_module += coeff* TensorModule.ext_multiply(FA(*word).coproduct()*unit_tensor_rc_module, graph_module)
    #         pickle_module += coeff* TensorModule.ext_multiply(FA(*word).change_basis(SchubertBasis), graph_module)

    #ring = (FA@FA)@Sx
    ring = (FA@FA)@(ASx@Sx)
    ring2 = (ASx@ASx)@Sx
    stink = ASx@Sx
    solution_element = ring.zero
    solution_element2 = ring2.zero
    solution_pickle = stink.zero

    def ring_elem_word(rc):
        return (ASx@Sx).ext_multiply(ASx(rc.perm,len(rc)), Sx(rc.polyvalue(x)))

    def ring_elem(rc):
        return (ASx@Sx).ext_multiply(ASx(rc.perm,len(rc)), Sx(rc.perm))

    solution_module = TensorModule()
    # THIS IS THE CORRECT COPRODUCT
    for seq in aseqs:
        #mod  = ASx(uncode(seq),n-1) * unit_rc_module
        perm = uncode(seq)
        elem = ASx(perm,n-1).change_basis(WordBasis)
        
        # accum_mod = TensorModule()
        # accum_elem = 0
        for word, coeff in elem.items():
            modmod = FA(*word).coproduct()*unit_tensor_rc_module
            for word_double, coeff_double in elem.items():
                mod = FA(*word_double) * unit_rc_module
                for rc, coeff2 in mod.items():
                    solution_module += coeff * coeff_double * coeff2 * TensorModule.ext_multiply((ASx@Sx)(((perm,n-1),rc.perm)), modmod)
            # accum_mod = TensorModule()
        
            #coeff * coeff2 * TensorModule.ext_multiply()

        # for perm in perms:
        #     perm_coeff  = ASx(perm,n-1).change_basis(WordBasis).get(seq,0)
        #     if perm_coeff != 0:
                # THIS IS IT
        # mod = FA(*seq) * unit_rc_module
        # for (rc1, rc2), coeff in mod.items():
        #     if len(rc1.perm) > n or len(rc2.perm) > n:
        #         continue
        #     solution_module2 += coeff * TensorModule.ext_multiply(TensorModule({RCGraphTensor(rc1, rc2): 1}), ASx)
        #     #
        # solution_module2 += TensorModule.ext_multiply(X,FA(*seq).coproduct()*unit_tensor_rc_module)
        #TensorModule.ext_multiply(FA(*word).coproduct()*unit_tensor_rc_module, FA(*word).change_basis(SchubertBasis)*unit_rc_module)

    # rabies = 0
    # for ((rc1, rc2), perm), coeff in solution_module.items():
    #     if len(rc1.perm) > n or len(rc2.perm) > n:
    #         continue
    #     rabies += coeff *((ASx@ASx)@Sx).ext_multiply((ASx@ASx).ext_multiply(ASx(rc1.perm,n-1), ASx(rc2.perm,n-1)),Sx(perm))
    # print(rabies)

    # print("TEST")

    rabies = 0
    R = (ASx@Sx)
    coprods = {}
    for (((perm, _),perm0), (rc1, rc2)), coeff in solution_module.items():
        #perm = rc0.perm
        if perm != perm0:
            print(f"{perm.trimcode} != {perm0.trimcode}")
            print(TensorModule({RCGraphTensor(rc1, rc2): coeff}))
            #continue
        #assert perm == perm0 or coeff == 0
        if len(rc1.perm) > n or len(rc2.perm) > n or len(perm) > n:
            continue
        # if product_too_big(rc1.perm, rc2.perm, n):
        #     continue
        #coprods[perm] = coprods.get(perm, 0) + coeff *(ASx@ASx).ext_multiply(ASx(rc1.perm, len(rc1)), ASx(rc2.perm, len(rc2)))
        coprods[perm0] = coprods.get(perm, 0) + coeff *(ASx@ASx).ext_multiply(ASx(rc1.perm, len(rc1)), ASx(rc2.perm, len(rc2)))
    
    for perm, elem in coprods.items():
        print(f"{perm.trimcode}")
        print(elem)
        check = ASx(perm, n-1).coproduct()
        print(check)
        diff = elem - check
        print(diff)
        assert all(v == 0 for v in diff.values()), f"Failed check on {perm}, diff = {diff}"
    exit()

    R = (ASx@Sx)
    for((rc1, rc2),rc), coeff in solution_module.items():
        solution_element = coeff*((R@R)@R).ext_multiply((R@R).ext_multiply(ring_elem(rc1, n), ring_elem(rc2, n)), ring_elem(rc,n))
        
        # if rc.is_principal:
        #     solution_element2 += coeff *ring2.ext_multiply((ASx@ASx).ext_multiply(FA(*seq1).change_basis(SchubertBasis),
        #                                                                         FA(*seq2).change_basis(SchubertBasis)), 
        #                                                                         Sx(rc.perm))
    for (schub,rc), coeff in pickle_module.items():
        solution_pickle += coeff * stink((schub, rc.perm))

    print(solution_element)
    # print(solution_pickle)
    # print(solution_element2)
    exit()
    # for seq in aseqs:
    #     pish_mod = FA(*seq) * unit_rc_module
    #     for rc, bob in pish_mod.items():
    #         solution_module3 += bob * TensorModule.ext_multiply(pish_mod,FA(*seq).coproduct()*unit_tensor_rc_module)
    #         if len(rc.perm) > n:
    #             continue
    #         #solution_module += coeff*TensorModule.ext_multiply(pish_mod,ASx(rc.perm, n-1).coproduct())
    #         # fist = ASx(rc.perm, n-1).coproduct()*unit_tensor_rc_module
    #         # for (rc1, rc2), coeff2 in fist.items():
    #         #     if vector_sum(rc1.length_vector(), rc2.length_vector()) == seq:
    #         #         #print(f"Adding {(rc, (rc1, rc2))} with coeff {bob*coeff*coeff2}")
    #         #solution_module2 += bob * TensorModule.ext_multiply(pish_mod,ASx(rc.perm,n-1).coproduct()*unit_tensor_rc_module)
    #         solution_module2 += bob * TensorModule.ext_multiply(pish_mod,ASx(rc.perm,n-1).coproduct()*unit_tensor_rc_module)
                    #solution_module2 += bob * coeff*TensorModule.ext_multiply(pish_mod,)
                #solution_module2 += coeff*TensorModule.ext_multiply(pish_mod,FA(*rc.length_vector()).coproduct()*unit_tensor_rc_module)

    # for perm in perms:
    #     pish_mod = (ASx(perm, n-1) * unit_rc_module)
    #     for rc, coeff in pish_mod.items():
    #         if len(rc.perm) > n:
    #             continue
    #         #solution_module += coeff*TensorModule.ext_multiply(pish_mod,ASx(rc.perm, n-1).coproduct())
    #         solution_module2 += coeff*TensorModule.ext_multiply(pish_mod,FA(*rc.length_vector()).coproduct()*unit_tensor_rc_module)

            # BIJECTION
    # print(solution_module2)
    # exit()
    tring = ASx@ASx
    # for seq in aseqs:
    #     schubelem1 = FA(*seq).change_basis(SchubertBasis)
    #     mod = RCGraphModule()
    #     for key, coeff in schubelem1.items():
    #         mod += coeff * ASx(*key) * unit_rc_module
    #     print(f"Sequence {seq}")
    #     print(mod)
    #     print(mod.schubvalue(Sx))
    # baby = 0
    # for perm in perms:
    #     baby += ASx(perm, n-1)
    
    # babyword = baby.change_basis(WordBasis)

    # for seq, coeff in babyword.items():
    #     pish_mod = (FA(*seq) * unit_rc_module)
    #     for rc, coeff2 in pish_mod.items():
    #         solution_module2 += coeff * coeff2 * TensorModule.ext_multiply(pish_mod, FreeAlgebraBasis.change_tensor_basis(ASx(rc.perm, n-1).coproduct(),WordBasis,WordBasis)*unit_tensor_rc_module)

        # cprd = schubelem1.coproduct()
        # for key0, coeff in schubelem1.items():
        #     for (key1, key2), coeff2 in cprd.items():
        #         solution_module2 += ASx(*key0).change_basis(WordBasis).get(seq,0)* coeff2 * TensorModule.ext_multiply(ASx(*key0)*unit_rc_module,tring((key1,key2))*unit_tensor_rc_module)

    #exit()
    products = {}
    #coproducts = {}
    coproducts2 = {}
    #coproducts3 = {}
    
    # for (rc0, (rc1, rc2)), coeff1 in solution_module.items():
    #     if len(rc0.perm) > n or len(rc1.perm) > n or len(rc2.perm) > n:
    #         continue
    #     if rc0.length_vector() == vector_sum(rc1.length_vector(), rc2.length_vector()):
    #         coproducts[rc0.perm] = coproducts.get(rc0.perm, 0) + coeff1 * tring((rc1.perm, rc2.perm))
    # for (rc0, (key1, key2)), coeff1 in solution_module.items():
    #     if len(rc0.perm) > n or len(key1[0]) > n or len(key2[0]) > n:
    #         continue
    #     coproducts[rc0.perm] = coproducts.get(rc0.perm, 0) + coeff1 * tring((key1, key2))

    rr0 = Sx@Sx
    rr = Sx@rr0
    sxt = rr.zero

    # for (rc0, (rc1, rc2)), coeff1 in solution_module2.items():
    #     if len(rc0.perm) > n or len(rc1.perm) > n or len(rc2.perm) > n:
    #         continue

        #if rc1.is_principal and rc2.is_principal:
        #coproducts2[rc0.perm] = coproducts2.get(rc0.perm, 0) + coeff1 * tring(((rc1.perm,n-1), (rc2.perm,n-1)))
        #sxt += coeff1*rr.ext_multiply(Sx(rc0.perm),rr0.ext_multiply(Sx(rc1.perm),Sx(rc2.perm)))

    

    for (rc0, (rc1, rc2)), coeff1 in solution_module2.items():
        if len(rc0.perm) > n or len(rc1.perm) > n or len(rc2.perm) > n:
            continue
        #coproducts2[rc0.perm] = coproducts2.get(rc0.perm, 0) + coeff1 * tring(((rc1.perm,n-1), (rc2.perm,n-1)))
        
        #shoeff = ASx(rc0.perm,n-1).change_basis(SchubertBasis).get(rc0.length_vector(),0)
        products[(rc1.perm, rc2.perm)] = products.get((rc1.perm, rc2.perm), 0) + coeff1 *  Sx(rc0.polyvalue(x))
        #if rc1.is_principal and rc2.is_principal:
            #products[(rc1.perm, rc2.perm)] = products.get((rc1.perm, rc2.perm), RCGraphModule()) + coeff1 * rc0
            #

    # for (perm1, perm2),elem in products.items():
    #     if any(len(perm) > n for perm in (Sx(perm1) * Sx(perm2)).keys()):
    #        continue 
    #         # check
    #     print(perm1.trimcode, perm2.trimcode)
    #     print(elem)
        #assert Sx(perm1) * Sx(perm2) == elem, f"Failed check {perm1}, {perm2} {elem=} {Sx(perm1) * Sx(perm2)=}"
    # for (rc0, (rc1, rc2)), coeff1 in solution_module3.items():
    #     if len(rc0.perm) > n or len(rc1.perm) > n or len(rc2.perm) > n:
    #         continue
    #     coproducts3[rc0.perm] = coproducts3.get(rc0.perm, 0) + coeff1 * tring(((rc1.perm,n-1), (rc2.perm,n-1)))
        # products[(key1[0], key2[0])] = products.get((key1[0], key2[0]), 0) + coeff1 * Sx(rc0.polyvalue(x))
        
    

    #productS = Sx([]).ring.zero
    for (perm1, (perm2, perm3)), coeff in sxt.items():
    # # rr = Sx@@x
    # # sxt = rr.zero
    #     products[(perm2, perm3)] = products.get((perm2, perm3), 0) + coeff * Sx(perm1)
        coproducts2[perm3] = coproducts2.get(perm3, 0) + coeff * (ASx@ASx).ext_multiply(ASx(perm2,n-1),ASx(perm2,n-1))
    
    for (perm1, perm2),elem in products.items():
        if any(len(perm) > n for perm in (Sx(perm1) * Sx(perm2)).keys()):
           continue 
        print(perm1.trimcode, perm2.trimcode)
        print(elem)

    #         # check
    #     assert Sx(perm1) * Sx(perm2) == elem, f"Failed check {perm1}, {perm2} {elem=} {Sx(perm1) * Sx(perm2)=}"

    for perm, elem in coproducts2.items():
        print(f"{perm.trimcode}")
        print(elem)
        print(ASx(perm, n-1).coproduct())


    # for perm, elem in coproducts2.items():
    #     print(f"{perm.trimcode}")
    #     print(elem)
    #     print(ASx(perm, n-1).coproduct())

    # for perm, elem in coproducts3.items():
    #     print(f"{perm.trimcode}")
    #     print(elem)
    #     print(ASx(perm, n-1).coproduct())

    exit()


if __name__ == "__main__":
    main()

