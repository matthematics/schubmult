import math
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

def len_coeff(rc):
    return FA(*rc.length_vector()).change_basis(SchubertBasis).get((rc.perm, len(rc)), 0)


def perm_coeff(perm, seq):
    return ASx(perm, len(seq)).change_basis(WordBasis).get(seq, 0)

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

def schub_inverse_pair(seq):
    schubseq = FA(*seq).change_basis(SchubertBasis)
    seqinv = Sx(expand_seq(seq,x))
    
    return schubseq, seqinv

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
    #ASx([]).ring @ Sx([]).ring
    test_addup = TensorModule()
    if n == 3:
        assert (2,0) in seqs
    result_dict = {}
    more_addup = TensorModule()
    aseqs = artin_sequences(n-1)
    assert len(aseqs) == math.factorial(n)
    upmod = TensorModule()
    for seq in aseqs:
                # left_result = TensorModule()
            #     addup = TensorModule()
                #the_elem = fa_elem_rc(seq, n)
        coprod = FA(*seq).coproduct()
        schubseq, seqinv = schub_inverse_pair(seq)
        for (seq1, seq2), coeff in coprod.items():
            schubseq1, seqinv1 = schub_inverse_pair(seq1)
            schubseq2, seqinv2 = schub_inverse_pair(seq2)
            upmod += coeff * TensorModule.ext_multiply(TensorModule.ext_multiply(schubseq,
                                                                         seqinv),
                                                    TensorModule.ext_multiply(
                                            TensorModule.ext_multiply(schubseq1, seqinv1),
                                            TensorModule.ext_multiply(schubseq2, seqinv2)))
            # upmod += coeff * TensorModule.ext_multiply(fa_elem_rc(seq, n), TensorModule.ext_multiply(fa_elem_rc(seq1, n), fa_elem_rc(seq2, n)))
            #upmod += coeff * TensorModule.ext_multiply(farp_elem_rc(seq, n), TensorModule.ext_multiply(farp_elem_rc(seq1, n), farp_elem_rc(seq2, n)))
            # rc after
            # upmod += coeff * TensorModule.ext_multiply(fa_elem_rc(seq, n), TensorModule.ext_multiply(fa_elem_rc(seq1, n), fa_elem_rc(seq2, n)))
    
    # fiddlemod = TensorModule()

    # for (((perm, _), seq), (((perm1, _), seq1), ((perm2, _), seq2))), coeff in upmod.items():
    #     # product = Sx(perm1) * Sx(perm2)
    #     # if any(len(permperm)>n for permperm in product.keys()):
    #     #     continue
    #     # # result_dict[(perm1, perm2)] = result_dict.get((perm1, perm2), RCGraphModule()) + coeff * rc
    #     # fiddlemod += coeff * TensorModule.ext_multiply(TensorModule.ext_multiply(Sx(perm),Sx(rc.polyvalue(x))),
    #     #                                                TensorModule.ext_multiply(TensorModule.ext_multiply(Sx(perm1), Sx(rc1.polyvalue(x))),
    #     #                                                                          TensorModule.ext_multiply(Sx(perm2), Sx(rc2.polyvalue(x)))))
    #     # seq_module = FA(*seq) * unit_rc_module
    #     # seq1_module = FA(*seq1) * unit_rc_module
    #     # seq2_module = FA(*seq2) * unit_rc_module
    #     # fiddlemod += coeff * TensorModule.ext_multiply(TensorModule.ext_multiply(Sx(perm),
    #     #                                                                          Sx([]).ring.from_dict({rc.perm: cff*perm_coeff(rc.perm, seq) for rc, cff in seq_module.items()})),
    #     #                                     TensorModule.ext_multiply(TensorModule.ext_multiply(Sx(perm1), 
    #     #                                                                                         Sx([]).ring.from_dict({rc.perm: cff*perm_coeff(rc.perm, seq1) for rc, cff in seq1_module.items()})),
    #     #                                                                 TensorModule.ext_multiply(Sx(perm2), 
    #     #                                                                                           Sx([]).ring.from_dict({rc.perm: cff*perm_coeff(rc.perm, seq2) for rc, cff in seq2_module.items()}))))
    #     fiddlemod += coeff * TensorModule.ext_multiply(TensorModule.ext_multiply(Sx(perm),
    #                                                                              Sx(expand_seq(seq, x))),
    #                                         TensorModule.ext_multiply(TensorModule.ext_multiply(Sx(perm1), 
    #                                                                                     Sx(expand_seq(seq1, x))),
    #                                                                     TensorModule.ext_multiply(Sx(perm2), 
    #                                                                                 Sx(expand_seq(seq2, x)))))
            
        #assert product.get(perm, 0) == coeff, f"{coeff=} {perm1.trimcode=} {perm2.trimcode=} {perm.trimcode=} {product=}"
    
    addups = {}
    for (((perm0, _), perm00), (((perm1, _), perm11), ((perm2, _), perm22))), coeff in upmod.items():
        
        

        if perm0 != perm00:
            continue
        if perm1 != perm11 or perm2 != perm22:
            continue
        addups[perm0] = addups.get(perm0, 0) + coeff * (ASx@ASx)(((perm1,n-1), (perm2,n-1)))
        # assert perm0 == rc.perm
        # assert perm1 == rc1.perm
        # assert perm2 == rc2.perm
        
        # if not rc1.is_principal and not rc2.is_principal:
        #     continue
        # assert coeff >= 0
        product = Sx(perm1) * Sx(perm2)
        # if rc.is_principal:
        #     if any(len(permperm)>n for permperm in product.keys()):
        #         continue
        assert product.get(perm0, 0) == coeff, f"{perm0.trimcode=} {perm0.trimcode=} {perm1.trimcode=} {perm2.trimcode=} {product=} {coeff=}"
        #     print(f"Success {perm0.trimcode=} {perm1.trimcode=} {perm2.trimcode=} {product=} {coeff=}")
    # for (perm1, perm2), mod in result_dict.items():
    #     product = Sx(perm1) * Sx(perm2)
    #     print(mod)
    #     assert product == Sx(mod.polyvalue(x))

    
    # cprod = {}

    # siphon_off = {}
    # tring = ASx([]).ring @ ASx([]).ring
    # da_module = TensorModule()
    # ring = (ASx([]).ring @ Sx([]).ring) @ (ASx([]).ring @ ASx([]).ring)
    # tester_module = ring.zero
    # for ((seq1, seq2), ((perm1, _), (perm2, _))), coeff in upmod.items():
    #     seq = vector_sum(seq1, seq2)
    #     elem = TensorModule.ext_multiply(FA(*seq).change_basis(SchubertBasis), Sx(expand_seq(seq, x)))
    #     # for rc, v in elem.items():
    #     #     siphon_off[(perm1, perm2)] = siphon_off.get((perm1, perm2), RCGraphModule()) + v * coeff * rc
    #     #da_module += coeff * TensorModule.ext_multiply(elem, tring(((perm1, n-1), (perm2, n-1))))
    #     da_module += coeff * TensorModule.ext_multiply(TensorModule.ext_multiply(elem,{(seq1, seq2): 1}), tring(((perm1, n-1), (perm2, n-1))))
    #     # tester_module += coeff * ring.ext_multiply(ring.rings[0].ext_multiply(FA(*seq1).change_basis(SchubertBasis),Sx(expand_seq(seq1,x))),
    #     #                                           ring.rings[1].ext_multiply(ASx(perm1,n-1), ASx(perm2,n-1)))

    # for (((perm0, _), perm00), ((perm1, _), (perm2, _))), coeff in tester_module.items():
    #     if coeff == 0:
    #         continue
    #     try:
    #         assert perm0 == perm00
    #     except AssertionError:
    #         continue
    #     product = Sx(perm1) * Sx(perm2)
    #     if any(len(permperm)>n for permperm in (product).keys()):
    #         continue
    #     assert product.get(perm00, 0) == coeff, f"Failure for {perm0.trimcode,perm1.trimcode,perm2.trimcode}: {product} at {perm00.trimcode} != {coeff}"

    #for (((perm, _), rc), ((perm1, _), (perm2, _))), coeff in da_module.items():
    # for ((((perm, _), perm0),(seq1, seq2)), ((perm1, _), (perm2, _))), coeff in da_module.items():
    #     if coeff == 0:
    #         continue
    #     #assert perm == perm0, f"{perm.trimcode} {perm0.trimcode} {perm1.trimcode} {perm2.trimcode} {coeff}"
    #     #siphon_off[(perm1, perm2)] = siphon_off.get((perm1, perm2), RCGraphModule()) + coeff * rc
    #     siphon_off[(perm1, perm2)] = siphon_off.get((perm1, perm2), 0) + coeff * Sx(perm0)
    #for (rc, (((perm1, _), rc1), ((perm2, _), rc2))), coeff in more_addup.items():
    # for perm, ring_elem in siphon_off.items():
    #     coproduct = ASx(perm, n-1).coproduct()
    #     #print(f"{perm.trimcode}: {ring_elem} vs {product}")
    #     assert all(v == 0 for v in (ring_elem - coproduct).values()), f"Failure for {perm.trimcode}: {ring_elem} != {coproduct}"
    #     print(f"Success for {perm.trimcode}: Sx({perm.trimcode})*Sx({perm.trimcode})")
    #     print(perm)
    # for (((perm, _), rc), ((perm1, _),(perm2, _))), coeff in upmod.items():
    #     #more_addup2[(rc1.perm, rc2.perm)] = more_addup2.get((rc1.perm, rc2.perm), RCGraphModule()) + coeff * rc
    #     #more_addup3[perm] = more_addup3.get(perm, 0) + coeff * (Sx([]).ring@Sx([]).ring)((perm1, perm2))
    #     assert rc.perm == perm
    #     cprod[(perm1, perm2)] = cprod.get((perm1, perm2), RCGraphModule()) + coeff * rc
    #     product = Sx(perm1) * Sx(perm2)
    #     #assert product.get(perm, 0) == coeff, f"Failure for {perm.trimcode}: {product.get(perm,0)} != {coeff}"
    #     #print(f"{perm.trimcode}: ({perm1.trimcode},{perm2.trimcode}) -> {coeff} Success")

    siphon_off = addups
    #for (perm1, perm2), coeff in siphon_off.items():
    for perm, coeff in siphon_off.items():
        #coproduct = Sx(perm1) * Sx(perm2)
        coproduct = ASx(perm, n-1).coproduct()
        try:
            assert all(v==S.Zero for v in (coproduct - coeff).values())
        except AssertionError:
            print(f"Failure for {perm.trimcode}: {coeff=} {coproduct=}")
            continue
        print(f"Success {perm.trimcode}")
    # for (perm1, perm2), mod in more_addup2.items():
    #     product = Sx(perm1) * Sx(perm2)
    #     if any(len(permperm)>n for permperm in (product).keys()):
    #         continue
    #     print(f"{perm1.trimcode,perm2.trimcode}")
    #     print(mod)
        # for rc in mod:
        #     assert product.get(rc.perm, 0) == coeff
    

        
    
    print("HELLS YEAH")
        #assert Sx(mod.polyvalue(x)) == product, f"Failure for {(perm1.trimcode, perm2.trimcode)}: {mod} != {product}"
        #assert mod == product, f"Failure for {(perm1.trimcode, perm2.trimcode)}: {mod} != {product}"
        #     module1 = FA(*seq).change_basis(SchubertBasis)
        #     #     #rc_graph_module_term = RCGraphModule({rc: v for rc, v in rc_graph_module_term.items() if len(rc.perm) <= n})
        #         #tmodule = TensorModule.ext_multiply(rc_graph_module_term, rc_graph_module_term)
        #     #         #if rc1 == rc2:
        #     # #                addup += TensorModule.ext_multiply(ASx(rc1.perm,len(rc1)), RCGraphModule({rc2: coeff}))

        #     coprod = FA(*seq).coproduct()
        #     #     #test_addup += addup
        #     #     #left_result = TensorModule.ext_multiply(addup, FreeAlgebraBasis.change_tensor_basis(coprod, SchubertBasis, SchubertBasis))
        #     #     left_result = addup
        #     #     # # zr_elem = zring(seq)
        #     #     # # xr_elem = xring(seq)
        #     #     # # module1 = TensorModule.ext_multiply(rc_graph_module_term, zr_elem)
        #     #     # # module2 = TensorModule.ext_multiply(xr_elem, rc_graph_module_term)
        #     #     # module1 = TensorModule.ext_multiply(rc_graph_module_term, rc_graph_module_term)
        #     #     # # module2 = TensorModule.ext_multiply(rc_graph_module_term, rc_graph_module_term)
        #     #     right_result = TensorModule()
        #     #     for (rc1, rc2), coeff00 in tmodule.items():
        #     #         for (seq1, seq2), val in coprod.items():
        #         #tmodule1 = seq_elem_rc(seq, n)
        #     #tmodule1 = schub_elem_rc(seq, n)
        #     tmodule1 = Sx(expand_seq(seq,x))
        #     for (seq1, seq2), val in coprod.items():
        #         #tmodule2 = fa_elem_rc(seq1, n)
        #         tmodule2 = schub_elem_rc(seq1, n)
        #         tmodule3 = schub_elem_rc(seq2, n)

        #         result += val * TensorModule.ext_multiply(tmodule1, TensorModule.ext_multiply(tmodule2, tmodule3))
        # print(result)
        #         #             schub_elem1 = FA(*seq1).change_basis(SchubertBasis)
        #         #             schub_elem2 = FA(*seq2).change_basis(SchubertBasis)
        #         #             #module2 = RCGraphModule({rc: v for rc, v in module2.items() if len(rc.perm) <= n})
        #         #             #module3 = RCGraphModule({rc: v for rc, v in module3.items() if len(rc.perm) <= n})
        #             #tmodule2 = TensorModule.ext_multiply(module2, module2)
        #             #tmodule3 = TensorModule.ext_multiply(module3, module3)

        #             # for rc1_left, coeff1 in module1.items():
        #             #     for rc2_left , coeff2 in module2.items():
        #             #         for rc3_left,  coeff3 in module3.items():
        #             #             new_coeff = coeff2 * coeff3 * coeff1*val
        #             #                 #result += new_coeff * TensorModule.ext_multiply(ASx(rc1_left.perm, len(rc1_left)),TensorModule({RCGraphTensor(rc2_left, rc3_left): 1}))
        #             #             result += new_coeff * (fa_elem_rc(rc1_left) @ (fa_elem_rc(rc2_left) @ fa_elem_rc(rc3_left)))
        #     #modules = {}
        #     #print(test_addup)
        #     # for ((perm, _), rc), coeff in test_addup.items():
        #     #     modules[perm] = modules.get(perm, RCGraphModule()) + RCGraphModule({rc: coeff})

        #     # for perm in modules:
        #     #     print(perm.trimcode)
        #     #     print(modules[perm])
        #     #exit()
        #     #             full_tmodule = TensorModule.ext_multiply(tmodule2,tmodule3)
        #     #             addup = TensorModule()
        #     #             for ((rc1_left, rc2_left),(rc1_right,rc2_right)), coeff0 in full_tmodule.items():
        #     #                 coeff = coeff0 * coeff00
        #     #                 if rc1 == rc2 and rc1_left == rc2_left and rc1_right == rc2_right:
        #     #                     #print(f'{rc1_left.perm.trimcode, rc1_right.perm.trimcode}, coeff={coeff}, val={val}, coeff1={coeff1}, coeff2={coeff2}, seq1={seq1}, seq2={seq2}')
        #     #                     # print(f"{schub_elem1=}, {schub_elem2=}")
        #     #                     new_coeff = coeff * schub_elem1.get((rc1_left.perm, len(rc1_left)), 0) * schub_elem2.get((rc1_right.perm, len(rc1_right)), 0)
        #     #                     #print(f"{new_coeff=}")
        #     #                     #addup += TensorModule.ext_multiply(TensorModule.ext_multiply(RCGraphModule({rc1_left: new_coeff}), RCGraphModule({rc1_right: 1})), 
        #     #                     #                                    TensorModule.ext_multiply(ASx(rc1_left.perm,len(rc1_left)), ASx(rc1_right.perm,len(rc1_right))))
        #     #                     result += TensorModule.ext_multiply(TensorModule.ext_multiply(ASx(rc1.perm,len(rc1)), RCGraphModule({rc2: new_coeff})),
        #     #                                                         TensorModule.ext_multiply(TensorModule.ext_multiply(RCGraphModule({rc1_left: coeff}), RCGraphModule({rc1_right: 1})),
        #     #                                                         TensorModule.ext_multiply(ASx(rc1_left.perm,len(rc1_left)), ASx(rc1_right.perm,len(rc1_right)))))
        #                 #right_result += addup
        #         #                                                                                             TensorModule.ext_multiply(module3,module3)),
        #         #                            FreeAlgebraBasis.change_tensor_basis((FA @ FA)((seq1,seq2)),SchubertBasis,SchubertBasis)))
        #         #result += TensorModule.ext_multiply(left_result,right_result)
        #     # buildup_mod = {}
        #     # for ((perm, _), rc), coeff in test_addup.items():
        #     #     assert Sx(rc.polyvalue(x)) == Sx(perm), f"Failure for {perm.trimcode}"
        #     #     assert all(rc.perm == perm for rc in test_addup[rc].keys()), f"Failure for {perm.trimcode}"
        #     # left_graph == right_graph
        #     #result2 = result

        #     # result_dict = {}
        #     # for ((perm, _),rc), coeff in addup.items():
        #     #     result_dict[perm] = result_dict.get(perm, RCGraphModule()) + RCGraphModule({rc: coeff})

        #     # for perm in result_dict:
        #     #     assert Sx(result_dict[perm].polyvalue(x)) == Sx(perm), f"Failure for {perm.trimcode}"
        #     #     assert all(rc.perm == perm for rc in result_dict[perm].keys()), f"Failure for {perm.trimcode}"
        #     # exit()
        #     # new_result = TensorModule()
        # failed = False
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
        
        # #for (((perm, _), seq), (((left_coprod_perm, _), seq1), ((right_coprod_perm, _), seq2))), value in result.items():
        # #for (((perm, _), rc_result), (((left_coprod_perm, _), rc1), ((right_coprod_perm, _), rc2))), value in result.items():
        # for (rc_result_right, ((rc1_left, rc1_right), (rc2_left, rc2_right))), value in result.items():
        #     left_coprod_perm, right_coprod_perm = rc1_right, rc2_right

        #     if any(len(permperm) > n for permperm, val in (Sx(left_coprod_perm) * Sx(right_coprod_perm)).items() if val != S.Zero):
        #         continue
            
        #     perm_pairs_total.add((left_coprod_perm, right_coprod_perm))
            
        #     # left_coprod_perm, right_coprod_perm = rc1.perm, rc2.perm
        #     # if len(perm) > n or len(left_coprod_perm) > n or len(right_coprod_perm) > n:
        #     #     continue
            
        #     #if rc_result_left.perm != rc_result_right.perm or rc1_left.perm != rc1_right.perm or rc2_left.perm != rc2_right.perm or not rc_result_left.is_principal or not rc1_right.is_principal or not rc2_right.is_principal:
        #     # assert rc_result_left.perm == rc_result_right.perm
        #     # assert rc1_left.perm == rc1_right.perm
        #     # assert rc2_left.perm == rc2_right.perm
        #     # if rc1_right != rc1_left:
        #     #     continue
        #     # if rc2_right != rc2_left:
        #     #     continue
        #     # if not rc1_right.is_principal:
        #     #     continue
        #     # if not rc2_right.is_principal:
        #     #     continue
        #     # if vector_sum(rc1_right.length_vector(), rc2_right.length_vector()) != rc_result_right.length_vector():
        #     #     continue
        #     # if not rc_result_right.is_principal:
        #     #     continue
        #     product = Sx(left_coprod_perm)*Sx(right_coprod_perm)
            
            

        #     # coeff = 1
        #     # # coeff = ASx(perm, len(rc_result)).change_basis(WordBasis).get(rc_result.length_vector(), 0) * ASx(left_coprod_perm, len(rc1)).change_basis(WordBasis).get(rc1.length_vector(), 0) * ASx(right_coprod_perm, len(rc2)).change_basis(WordBasis).get(rc2.length_vector(), 0)
        #     # # if coeff == 0:
        #     # #     continue
        #     # perm_pairs_total.add((left_coprod_perm, right_coprod_perm))
        #     # #addups[(left_coprod_perm, right_coprod_perm)] = addups.get((left_coprod_perm, right_coprod_perm), RCGraphModule()) + RCGraphModule({result_rc: value})
        #     # if rc1.perm != left_coprod_perm or rc2.perm != right_coprod_perm or rc_result.perm != perm:
        #     #     continue
        #     product = Sx(left_coprod_perm) * Sx(right_coprod_perm)
        #     perm = rc_result_right
        #     #RCGraphModule(dict.fromkeys(RCGraph.all_rc_graphs(result_rc.perm, n-1), value))
        #     #if perm != uncode(seq) or left_coprod_perm != uncode(seq1) or right_coprod_perm != uncode(seq2):
        #     # if perm != rc_result.perm or left_coprod_perm != rc1.perm or right_coprod_perm != rc2.perm:
        #     #     continue
        #     # perm = rc_result.perm
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
        #     # #assert perm == uncode(vector_sum(rc1.length_vector(), rc2.length_vector())), f"Failure for {(left_coprod_perm.trimcode, right_coprod_perm.trimcode)}"
        #     # #addups[(rc1, rc2)] = addups.get((rc1, rc2), 0) + value * Sx(perm)
        #     # #addups[(rc1, rc2)] = addups.get((rc1, rc2), RCGraphModule()) + value * rc_result
        #     addups[(left_coprod_perm, right_coprod_perm)] = addups.get((left_coprod_perm, right_coprod_perm), 0) + value * Sx(perm)

        # # failed = {}

        # # for (perm1, perm2), val in addups.items():
        # #     #perm1, perm2 = rc1.perm, rc2.perm
        # #     product = Sx(perm1) * Sx(perm2)
        # #     print("perm1, perm2 =", perm1.trimcode, perm2.trimcode)
        # #     print(val)
        # #     valu = Sx(val.polyvalue(x))
        # #     try:
        # #         assert all(v == 0 for v in (product - valu).values()), f"Failure for {(perm1.trimcode, perm2.trimcode)}: {val} != {Sx(perm1)*Sx(perm2)}"
        # #         print(f"Success for {(perm1.trimcode, perm2.trimcode)}: Sx({perm1.trimcode})*Sx({perm2.trimcode})")
        # #         failed[(perm1, perm2)] = False
        # #     except Exception as e:
        # #         print(f"Failure for {(perm1.trimcode, perm2.trimcode)}: {valu} != {Sx(perm1)*Sx(perm2)}")
        # #         if (perm1, perm2) not in failed:
        # #             failed[(perm1, perm2)] = True
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
        
        # for (perm1, perm2), module in addups.items():
        #     print(f"THIS IS {(perm1.trimcode, perm2.trimcode)}")
        #     print(module)
        #     assert module == Sx(perm1) * Sx(perm2)
    

            

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

        # success = 0
        # #if all(not v for v in failed.values()):
        # if not failed:
        #     print(f"YEAH!!! {len(perms)=} {len(perm_pairs)=} {len(perm_pairs_total)=}", file=sys.stderr)
        #     print(f"YEAH!!! {len(perms)=} {len(perm_pairs)=} {len(perm_pairs_total)=}", file=sys.stderr)
        # # else:
        # #     for (p1,p2), v in failed.items():
        # #         if v:
        # #             print(f"FAILED FOR {(p1.trimcode,p2.trimcode)}", file=sys.stderr)
        # #         else:
        # #             success += 1
        # #     print("Number of successes:", success, file=sys.stderr)

if __name__ == "__main__":
    main()
