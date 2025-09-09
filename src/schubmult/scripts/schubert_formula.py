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
    #ASx([]).ring @ Sx([]).ring
    test_addup = TensorModule()
    if n == 3:
        assert (0,1) in seqs
    result_dict = {}
    more_addup = TensorModule()
    aseqs = artin_sequences(n-1)
    upmod = TensorModule()
    for seq in aseqs:
        # upmod += TensorModule.ext_multiply(ng_elem_rc(seq,n),
        #                                    FreeAlgebraBasis.change_tensor_basis(FA(*seq).coproduct(),SchubertBasis,SchubertBasis))
        upmod += TensorModule.ext_multiply(TensorModule.ext_multiply(FA(*seq).change_basis(SchubertBasis),{seq: 1}),
                                           FreeAlgebraBasis.change_tensor_basis(FA(*seq).coproduct(),SchubertBasis,SchubertBasis))

    addup = {}
    # upmod2 = TensorModule()
    # for (seq, (perm, ((perm1, _), (perm2, _)))), coeff in upmod.items():
    #     upmod2 += coeff * TensorModule.ext_multiply(fa_elem_rc(seq,n),TensorModule.ext_multiply(Sx(perm), 
    #                                        (ASx@ASx)(((perm1, n-1),(perm2,n-1)))))

    #for (((rc11, rc1)), (((perm1, _), (perm2, _)))), coeff in upmod.items():
    for (((perm0, _), seq), (((perm1, _), (perm2, _)))), coeff in upmod.items():
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
        addup[(perm1, perm2)] = addup.get((perm1, perm2), 0) + coeff * perm_coeff(perm0, seq)*Sx(expand_seq(seq,x))

    for (perm1, perm2), coeff in addup.items():
        product = Sx(perm1) * Sx(perm2)
        if any(len(perm) > n for perm in product.keys()):
            continue
        try:
            assert product == coeff
            #assert product == Sx(coeff.polyvalue(x))
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
