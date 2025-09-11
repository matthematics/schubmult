import sys

from schubmult import ASx, Permutation, uncode
from schubmult.abc import x
from schubmult.rings import FA, Sx, WordBasis
from schubmult.rings.free_algebra_basis import FreeAlgebraBasis, SchubertBasis
from schubmult.rings.rc_graph_module import RCGraph, RCGraphModule, RCGraphTensor, TensorModule
from schubmult.utils.perm_utils import artin_sequences


def vector_sum(v1, v2):
    return tuple([a + b for a, b in zip(v1, v2)])


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
    for perm in perms:
        pish_mod = (ASx(perm, n-1) * unit_rc_module)
        for rc, coeff in pish_mod.items():
            if len(rc.perm) > n:
                continue
            #solution_module += coeff*TensorModule.ext_multiply(pish_mod,ASx(rc.perm, n-1).coproduct())
            solution_module2 += coeff*TensorModule.ext_multiply(pish_mod,ASx(rc.perm, n-1).coproduct()*unit_tensor_rc_module)

    # for perm in perms:
    #     pish_mod = (ASx(perm, n-1) * unit_rc_module)
    #     for rc, coeff in pish_mod.items():
    #         if len(rc.perm) > n:
    #             continue
    #         #solution_module += coeff*TensorModule.ext_multiply(pish_mod,ASx(rc.perm, n-1).coproduct())
    #         solution_module2 += coeff*TensorModule.ext_multiply(pish_mod,FA(*rc.length_vector()).coproduct()*unit_tensor_rc_module)

            # BIJECTION

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

    for (rc0, (rc1, rc2)), coeff1 in solution_module2.items():
        if len(rc0.perm) > n or len(rc1.perm) > n or len(rc2.perm) > n:
            continue
        coproducts2[rc0.perm] = coproducts2.get(rc0.perm, 0) + coeff1 * tring(((rc1.perm,n-1), (rc2.perm,n-1)))
        #products[(rc1.perm, rc2.perm)] = products.get((rc1.perm, rc2.perm), 0) + coeff1 * Sx(rc0.polyvalue(x))

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
        
    # for (perm1, perm2),elem in products.items():
    #     if any(len(perm) > n for perm in (Sx(perm1) * Sx(perm2)).keys()):
    #        continue 
    #         # check
    #     assert Sx(perm1) * Sx(perm2) == elem, f"Failed check {perm1}, {perm2} {elem=} {Sx(perm1) * Sx(perm2)=}"

    # for perm, elem in coproducts.items():
    #     print(f"{perm.trimcode}")
    #     print(elem)
    #     print(ASx(perm, n-1).coproduct())


    for perm, elem in coproducts2.items():
        print(f"{perm.trimcode}")
        print(elem)
        print(ASx(perm, n-1).coproduct())

    # for perm, elem in coproducts3.items():
    #     print(f"{perm.trimcode}")
    #     print(elem)
    #     print(ASx(perm, n-1).coproduct())

    exit()


if __name__ == "__main__":
    main()
