import sys

from schubmult import ASx, Permutation, uncode
from schubmult.abc import x, y, z
from schubmult.rings import FA, FreeAlgebra, FreeAlgebraBasis, MonomialBasis, NilHeckeRing, PolynomialAlgebra, SchubertBasis, SingleSchubertRing, Sx, TensorRing, WordBasis
from schubmult.rings.rc_graph_module import RCGraph, RCGraphModule, RCGraphTensor, TensorModule, all_fa_degree
from schubmult.symbolic import S, expand, expand_seq
from schubmult.utils.perm_utils import artin_sequences


def vector_sum(v1, v2):
    return tuple([a + b for a, b in zip(v1, v2)])

# def fa_elem_rc(rc):
#     #ag = PolynomialAlgebra(MonomialBasis(x, len(rc.perm)))
#     coeff = ASx(rc.perm, len(rc)).change_basis(WordBasis).get(rc.length_vector(), 0)
#     #return TensorModule.ext_multiply(RCGraphModule({rc: 1}), ag(*rc.length_vector()))
#     return TensorModule.ext_multiply(FA(*rc.length_vector()).change_basis(SchubertBasis), RCGraphModule({rc: coeff}))

# def rc_elem_rc(rc):
#     #ag = PolynomialAlgebra(MonomialBasis(x, len(rc.perm)))
#     coeff = ASx(rc.perm, len(rc)).change_basis(WordBasis).get(rc.length_vector(), 0)
#     #return TensorModule.ext_multiply(RCGraphModule({rc: 1}), ag(*rc.length_vector()))
#     return TensorModule.ext_multiply(RCGraphModule({rc: 1}), RCGraphModule({rc: coeff}))

def rc_elem_rc(seq, n):
    graphs = FA(*seq) * RCGraphModule({RCGraph(): 1})
    #ag = PolynomialAlgebra(MonomialBasis(x, len(rc.perm)))
    addup = TensorModule()
    for rc, coeff in graphs.items():
        if len(rc.perm) > n:
            return TensorModule()
        left_mod = RCGraphModule({RCGraph(): 1})
        left_mod = ASx(rc.perm, len(rc))# * left_mod
        #left_mod = RCGraphModule({r: v for r, v in left_mod.items() if len(r.perm) <= n})
        addup += TensorModule.ext_multiply(left_mod, RCGraphModule({rc: left_mod.change_basis(WordBasis).get(seq, 0) * coeff}))
    return addup

def seq_elem_rc(seq, n):
    #ag = PolynomialAlgebra(MonomialBasis(x, len(rc.perm)))
    mod1 = FA(*seq) * RCGraphModule({RCGraph(): 1})
    
    ret_mod = TensorModule()
    for rc, coeff in mod1.items():
        if len(rc.perm) > n:
            continue
        ret_mod += TensorModule.ext_multiply(RCGraphModule({rc: 1}),RCGraphModule(dict.fromkeys(RCGraph.all_rc_graphs(rc.perm, len(seq)), coeff)))

    # for (rc1, rc2), val in mod.items():
    #     if len(rc1.perm) > n or len(rc2.perm) > n:
    #         continue
    #     coeff = ASx(rc1.perm, len(rc1)).change_basis(WordBasis).get(rc2.length_vector(), 0)
    #     ret_mod += TensorModule.ext_multiply(RCGraphModule({rc1: val}), RCGraphModule({rc2: coeff}))
    return ret_mod

def len_coeff(perm, seq):
    return FA(*seq).change_basis(SchubertBasis).get((perm, len(seq)), 0)


def perm_coeff(perm, seq):
    return ASx(perm, len(seq)).change_basis(WordBasis).get(seq, 0)

def sx_coeff(perm, seq):
    return Sx(expand_seq(seq,x)).get(perm, 0)

def full_coeff(rc):
    return len_coeff(rc) * perm_coeff(rc)

def fa_elem_rc(seq, n):
    #ag = PolynomialAlgebra(MonomialBasis(x, len(rc.perm)))
    mod1 = FA(*seq) * RCGraphModule({RCGraph(): 1})
    
    ret_mod = TensorModule()
    for rc, coeff in mod1.items():
        if len(rc.perm) > n:
            continue
        ret_mod += TensorModule.ext_multiply(ASx(rc.perm, len(rc)),RCGraphModule({rc: coeff}))

    # for (rc1, rc2), val in mod.items():
    #     if len(rc1.perm) > n or len(rc2.perm) > n:
    #         continue
    #     coeff = ASx(rc1.perm, len(rc1)).change_basis(WordBasis).get(rc2.length_vector(), 0)
    #     ret_mod += TensorModule.ext_multiply(RCGraphModule({rc1: val}), RCGraphModule({rc2: coeff}))
    return ret_mod

def ng_elem_rc(seq, n):
    #ag = PolynomialAlgebra(MonomialBasis(x, len(rc.perm)))
    return TensorModule.ext_multiply(FA(*seq).change_basis(SchubertBasis), {seq: 1})

    # for (/rc1, rc2), val in mod.items():
    #     if len(rc1.perm) > n or len(rc2.perm) > n:
    #         continue
    #     coeff = ASx(rc1.perm, len(rc1)).change_basis(WordBasis).get(rc2.length_vector(), 0)
    #     ret_mod += TensorModule.ext_multiply(RCGraphModule({rc1: val}), RCGraphModule({rc2: coeff}))
    #return ret_mod


def farp_elem_rc(seq, n):
    #ag = PolynomialAlgebra(MonomialBasis(x, len(rc.perm)))
    mod1 = FA(*seq) * RCGraphModule({RCGraph(): 1})
    
    ret_mod = TensorModule()
    for rc, coeff in mod1.items():
        if len(rc.perm) > n:
            continue
        ret_mod += TensorModule.ext_multiply(ASx(rc.perm, len(rc)),RCGraphModule({rc: coeff*perm_coeff(rc)}))

    # for (rc1, rc2), val in mod.items():
    #     if len(rc1.perm) > n or len(rc2.perm) > n:
    #         continue
    #     coeff = ASx(rc1.perm, len(rc1)).change_basis(WordBasis).get(rc2.length_vector(), 0)
    #     ret_mod += TensorModule.ext_multiply(RCGraphModule({rc1: val}), RCGraphModule({rc2: coeff}))
    return ret_mod

def schub_elem_rc(seq, n):
    #ag = PolynomialAlgebra(MonomialBasis(x, len(rc.perm)))
    mod1 = FA(*seq) * RCGraphModule({RCGraph(): 1})
    
    ret_mod = TensorModule()
    for rc, coeff in mod1.items():
        if len(rc.perm) > n:
            continue
        #ret_mod += TensorModule.ext_multiply(ASx(rc.perm, len(rc)),Sx(RCGraphModule({rc: coeff}).polyvalue(x)))
        ret_mod += coeff * TensorModule.ext_multiply(ASx(rc.perm, len(rc)),Sx(expand_seq(rc.length_vector(),x)))

    # for (rc1, rc2), val in mod.items():
    #     if len(rc1.perm) > n or len(rc2.perm) > n:
    #         continue
    #     coeff = ASx(rc1.perm, len(rc1)).change_basis(WordBasis).get(rc2.length_vector(), 0)
    #     ret_mod += TensorModule.ext_multiply(RCGraphModule({rc1: val}), RCGraphModule({rc2: coeff}))
    return ret_mod


def main():
    n = int(sys.argv[1])

    seqs = artin_sequences(n - 1)
    xring = PolynomialAlgebra(MonomialBasis(x, n - 1))

    yring = PolynomialAlgebra(MonomialBasis(y,n-1))
    zring = PolynomialAlgebra(MonomialBasis(z,n-1))

    unit_rc_module = RCGraphModule({RCGraph(): 1})
    # result = TensorModule()
    result = 0
    degree = (n*(n-1))//2
    # 100% positive!
    seqs = set()
    for deg in range(degree+1):
        seqs.update(all_fa_degree(deg, n-1))
    #ASx([]).ring @ Sx([]).ring
    test_addup = TensorModule()
    if n == 3:
        assert (0,1) in seqs
    result_dict = {}
    more_addup = TensorModule()
    aseqs = artin_sequences(n-1)
    upmod = TensorModule()
    upmod2 = TensorModule()
    for seq in seqs:
        # upmod += TensorModule.ext_multiply(ng_elem_rc(seq,n),
        #                                    FreeAlgebraBasis.change_tensor_basis(FA(*seq).coproduct(),SchubertBasis,SchubertBasis))
        # upmod += TensorModule.ext_multiply(FA(*seq)*unit_rc_module,
        #                                   FreeAlgebraBasis.change_tensor_basis(FA(*seq).coproduct(),SchubertBasis,SchubertBasis))
        upmod2 += TensorModule.ext_multiply(FA(*seq).change_basis(SchubertBasis),
                                            FA(*seq).coproduct()*TensorModule.ext_multiply(unit_rc_module,unit_rc_module))
    addup = {}
    addup2 = {}
    addup0 = {}
    ring = ASx@ASx
    
    def rcs(seq, perm):
        return RCGraphModule({rc: v for rc, v in (FA(*seq)*RCGraphModule({RCGraph(): 1})).items() if rc.perm == perm})

    for (key, (rc1, rc2)), coeff in upmod2.items():
        if len(key[0]) > n:# or len(rc1.perm) > n or len(rc2.perm) > n:
            continue
        #addup[rc0.perm] = addup.get(rc0.perm, 0) + ASx(rc0.perm, n-1).change_basis(WordBasis).get(rc0.length_vector(),0)*coeff * ring((rc1, rc2))
        seq = vector_sum(rc1.length_vector(),rc2.length_vector())
        # USE THIS TO JUST ADD AN RC
        addup[(rc1.perm,rc2.perm)] = addup.get((rc1.perm,rc2.perm), 0) + ASx(*key).change_basis(WordBasis).get(seq,0)*coeff *Sx(expand_seq(seq,x)) 
        #addup[(rc1.perm,rc2.perm)] = addup.get((rc1.perm,rc2.perm), RCGraphModule()) + ASx(*key).change_basis(WordBasis).get(seq,0)*coeff *rcs(seq, key[0])
        #addup[key[0]] = addup.get(key[0], ring.zero) + perm_coeff(rc1.perm,rc1.length_vector())*perm_coeff(rc2.perm,rc2.length_vector())*coeff *ring(((rc1.perm,len(rc1)),(rc2.perm,len(rc2))))
        #Sx(expand_seq(seq,x)) 
        #ring(rc1.polyvalue(x))*Sx(rc2.polyvalue(x))
        #addup[(rc1[0],rc2[0])] = addup.get((rc1[0],rc2[0]), RCGraphModule()) + ASx(rc0.perm, n-1).change_basis(WordBasis).get(rc0.length_vector(),0)*coeff * rc0

    # for (rc0, (rc1, rc2)), coeff in upmod.items():
    #     if len(rc0.perm) > n or len(rc1[0]) > n or len(rc2[0]) > n:
    #         continue
    #     #addup[rc0.perm] = addup.get(rc0.perm, 0) + ASx(rc0.perm, n-1).change_basis(WordBasis).get(rc0.length_vector(),0)*coeff * ring((rc1, rc2))
    #     addup0[(rc1[0],rc2[0])] = addup0.get((rc1[0],rc2[0]), 0) + ASx(rc0.perm, n-1).change_basis(WordBasis).get(rc0.length_vector(),0)*coeff * Sx(rc0.polyvalue(x))
    #     addup[(rc1[0],rc2[0])] = addup.get((rc1[0],rc2[0]), RCGraphModule()) + ASx(rc0.perm, n-1).change_basis(WordBasis).get(rc0.length_vector(),0)*coeff * rc0
    #     #ASx(rc0.perm, n-1).change_basis(WordBasis).get(rc0.length_vector(),0)*coeff * ring((rc1, rc2))
    #     #Sx(expand_seq(rc0.length_vector(),x)).get(rc0.perm,0)*coeff * ring((rc1, rc2))

    #for perm, elem in addup.items():
    #     coproduct = ASx(perm, n-1).coproduct()
    #     if any(len(perm1[0][0]) > n or len(perm1[1][0])>n for perm1 in coproduct.keys() if coproduct[perm1]!=0):
    #         continue
    #     print(perm.trimcode)
    #     print(elem)
    #     try:
    #         #assert product == elem
    #         assert all(v == 0 for v in (coproduct-elem).values())
            
    #     except AssertionError:
    #         print(f"Failure {perm.trimcode}")
    #         print("Expected")
    #         print(coproduct)
    #         print("Got")
    #         print(elem)
    #         exit()
    #     print(f"Success {perm.trimcode}")
    #     print(elem)
    # print(len(addup))
    # exit()
    # #     assert product == addup0[(perm1, perm2)], f"Failure on {perm1.trimcode} {perm2.trimcode}\nExpected {product}\nGot {addup0[(perm1, perm2)]}"

    # # for (key, (rc1, rc2)), coeff in upmod2.items():
    # #     #print(key, rc1, rc2, coeff)
    # #     if len(rc1.perm) > n or len(rc2.perm) > n:
    # #         continue
    # #     if len(key[0]) > n:
    # #         continue
    # #     #perm_coeff(rc1.perm,rc1.length_vector())*perm_coeff(rc2.perm,rc2.length_vector())
    # #     seq = vector_sum(rc1.length_vector(),rc2.length_vector())
    # #     #addup2[(rc1.perm, rc2.perm)] = addup2.get((rc1.perm, rc2.perm), 0) + coeff * ASx(*key).change_basis(WordBasis).get(vector_sum(rc1.length_vector(),rc2.length_vector()),0)*Sx(key[0])
    # #     addup[(rc1.perm, rc2.perm)] = addup.get((rc1.perm, rc2.perm), 0) + coeff * ASx(*key).change_basis(WordBasis).get(seq,0)*Sx(key[0])
    # #     #RCGraphModule({k: v for k, v in (FA(*seq)*unit_rc_module).items() if k.perm == key[0]})
    # #     #ASx(*key).change_basis(WordBasis).get(vector_sum(rc1.length_vector(),rc2.length_vector()),0)*ring(((rc1.perm,len(rc1)),(rc2.perm,len(rc2))))
    # # # for perm, elem in addup2.items():
    # # #     #perm = rc0.perm
    # # #     #perm = rc0
    # # #     diff = ASx(perm, n-1).coproduct() - elem
    # # #     assert all(v == 0 for v in diff.values()), f"Failure on {perm.trimcode}\nExpected {ASx(perm, n-1).coproduct()}\nGot {elem}"
    for (perm1, perm2), elem in addup.items():
        product = Sx(perm1) * Sx(perm2)
        if any(len(perm) > n for perm in product.keys()):
            continue
        try:
            #assert product == elem
            assert product == addup[(perm1, perm2)]
            print(addup[(perm1, perm2)])
        except AssertionError:
            print(f"Failure {perm1.trimcode} {perm2.trimcode}")
            print("Expected")
            print(product)
            print("Got")
            print(elem)
            
        print(f"Success {perm1.trimcode} {perm2.trimcode}")
        
    # for perm, elem in addup2.items():
    #     coproduct = ASx(perm, n-1).coproduct()
    #     if any(len(perm1[0][0]) > n or len(perm1[1][0])>n for perm1 in coproduct.keys()):
    #         continue
    #     try:
    #         #assert product == elem
    #         assert all(v == 0 for k,v in (coproduct - elem).items())
    #     except AssertionError:
    #         print(f"Failure {perm.trimcode}")
    #         print("Expected")
    #         print(coproduct)
    #         print("Got")
    #         print(elem)
    #         continue
    #     print(f"Success {perm.trimcode}")
    #     print(elem)

    exit()
    # upmod2 = TensorModule()
    # for (seq, (perm, ((perm1, _), (perm2, _)))), coeff in upmod.items():
    #     upmod2 += coeff * TensorModule.ext_multiply(fa_elem_rc(seq,n),TensorModule.ext_multiply(Sx(perm), 
    #                                        (ASx@ASx)(((perm1, n-1),(perm2,n-1)))))

    #for (((rc11, rc1)), (((perm1, _), (perm2, _)))), coeff in upmod.items():
    #for (((perm0, _), rc1), (((perm1, _), (perm2, _)))), coeff in upmod.items():
    upmod2 = TensorModule()
    for (seq, (((perm1, _), (perm2, _)))), coeff in upmod.items():
        #addup[(perm1, perm2)] = addup.get((perm1, perm2), 0) + coeff * rc1
        # elem1 = ASx(perm1, n-1).change_basis(WordBasis)
        # elem2 = ASx(perm2, n-1).change_basis(WordBasis)
        # sum_elem = FA.zero
        # for seq1, coeff1 in elem1.items():
        #     for seq2, coeff2 in elem2.items():
        #         new_seq = vector_sum(seq1, seq2)
        #         sum_elem += coeff1 * coeff2 * FA(*new_seq)
        # print(sum_elem)
        #if sum_elem.get(rc1.length_vector(), 0) != 0:
        # addup[(perm1, perm2)] = addup.get((perm1, perm2), RCGraphModule()) + coeff *perm_coeff(perm0, rc1.length_vector())*rc1
        # addup[(perm1, perm2)] = addup2.get((perm1, perm2), RCGraphModule()) + coeff * rc1
        upmod2 += coeff * TensorModule.ext_multiply(fa_elem_rc(seq,n),(Sx@Sx)((perm1,perm2)))

    for (((perm0, _), rc1), ((perm1, perm2))), coeff in upmod2.items():
        addup[(perm1, perm2)] = addup.get((perm1, perm2), RCGraphModule()) + coeff * rc1

    # for (perm0, perm1, perm2), coeff in addup2.items():
    #     if all(rc in coeff.keys() for rc in RCGraph.all_rc_graphs(perm0, n-1)):
    #         addup[(perm1, perm2)] = addup.get((perm1, perm2), RCGraphModule()) + coeff

    for (perm1, perm2), coeff in addup.items():
        product = Sx(perm1) * Sx(perm2)
        if any(len(perm) > n for perm in product.keys()):
            continue
        try:
            #assert product == coeff
            assert product == Sx(coeff.polyvalue(x))
        except AssertionError:
            print(f"Failure {perm1.trimcode} {perm2.trimcode}")
            print("Expected")
            print(product)
            print("Got")
            print(coeff)
            continue
        print(f"Success {perm1.trimcode} {perm2.trimcode}")
        print(coeff)

if __name__ == "__main__":
    main()
