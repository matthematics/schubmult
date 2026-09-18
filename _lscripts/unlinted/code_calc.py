from schubmult import *
from schubmult.utils.tuple_utils import pad_tuple
from functools import cache
from schubmult.combinatorics.crystal_graph import CrystalGraphTensor
from schubmult.rings.combinatorial.hw_rc_ring import HWRCGraphRing
r = HWRCGraphRing()
CB = FreeAlgebra(CompositionSchubertBasis)

def crystal_coprod(perm, length):
    result = 0
    the_prin = RCGraph.principal_rc(perm, length)
    the_prin_hw, raise_seq = the_prin.to_highest_weight()
    # first_word = FA(*the_prin.length_vector)
    first_word = CB(*the_prin.length_vector).change_basis(WordBasis)
    #the_prin.to_free_algebra_element().coproduct()
    for word, coeff in first_word.items():
        word_coprod = FA(*word).coproduct()
        for (word1, word2), coeff2 in word_coprod.items():
            the_tensors = r.monomial(*word1) @ r.monomial(*word2)
            for (rc1, rc2), coeff3 in the_tensors.items():
                # tensor = CrystalGraphTensor(rc1, rc2)
                # if tensor.crystal_weight == the_prin_hw.length_vector:        
                result += coeff * coeff2 * coeff3 * r(rc1) @ r(rc2)
    #assert all(v >= 0 for _, v in result.items()), f"Negative coefficient in coproduct of {perm.trimcode}, got {result}"
    return result

@cache
def coproduct(perm, length, the_weight=None):
    from schubmult.combinatorics.crystal_graph import CrystalGraphTensor
    if perm.inv == 0:
        one =  r(RCGraph([()])).resize(length)
        return one @ one

    result = 0
    # if length == 2:
    #     trm_code = perm.trimcode
    #     if len(trm_code) == 1:
    #         trm_code = (trm_code[0], 0)
    #     result = 0
    #     for i in range(trm_code[0] + 1):
    #         for j in range(trm_code[1] + 1):
    #             rc_set = r.monomial(i, j) @ r.monomial(trm_code[0] - i, trm_code[1] - j)
    #             def decompose_tensor(rc1, rc2, coeff, flipped=False):
    #                 if flipped:
    #                     tensor = CrystalGraphTensor(rc2, rc1)
    #                 else:
    #                     tensor = CrystalGraphTensor(rc1, rc2)
    #                 # if perm.is_dominant:
    #                 #     if tensor.is_highest_weight and tensor.is_lowest_weight:
    #                 #         return coeff * ASx(rc1.perm, length) @ ASx(rc2.perm, length)
    #                 #     return 0
    #                 tensor_hw, _ = tensor.to_highest_weight()
    #                 if not perm.is_dominant and uncode(tuple(reversed(tensor_hw.crystal_weight))) == perm:
    #                     return coeff * ASx(rc1.perm, length) @ ASx(rc2.perm, length)
    #                 if tensor == tensor_hw and tensor_hw.is_lowest_weight and uncode(tensor.crystal_weight) == perm:
    #                     return coeff * ASx(rc1.perm, length) @ ASx(rc2.perm, length)
    #                 return r.zero
    #             for (rc1, rc2), coeff in rc_set.items():
    #                 #print(f"Decomposing tensor of {rc1} and {rc2} with coeff {coeff}")
    #                 if rc1.perm.inv == 0:
    #                     if rc2.perm == perm:
    #                         result += coeff * ASx(rc1.perm, length) @ ASx(rc2.perm, length)
    #                     continue
    #                 if rc2.perm.inv == 0:
    #                     if rc1.perm == perm:
    #                         result += coeff * ASx(rc1.perm, length) @ ASx(rc2.perm, length)
    #                     continue
    #                 if rc1.perm.is_dominant:
    #                     result += decompose_tensor(rc1, rc2, coeff)
    #                 elif rc2.perm.is_dominant:
    #                     # if rc1.perm.is_dominant:
    #                     #     if perm.is_dominant:
    #                     #         result += coeff * ASx(rc1.perm, length) @ ASx(rc2.perm, length)
    #                     #     continue
    #                     #flipped_rc1, flipped_rc2 = rc2, rc1
    #                     result += decompose_tensor(rc1, rc2, coeff, flipped=True)
    #                 else:
    #                     # if perm.is_dominant:
    #                     #     continue
    #                     result += decompose_tensor(rc1, rc2, coeff)
    #     return result
    #                 #result += coeff * r.from_rc_graph(rc).to_free_algebra_element() @ r.from_rc_graph(rc).to_free_algebra_element()
    #     #     for i in range(trm_code[0] + 1):
    #     #         result += r.monomial(i, 0).to_free_algebra_element() @ r.monomial(trm_code[0] - i, 0).to_free_algebra_element()
    #     #     result -= coproduct(uncode([0, trm_code[0]]), 2)
    #     #     return result
    #     # # result = r.monomial
    #     # if trm_code[0] == 0:
    #     #     for i in range(trm_code[1] + 1):
    #     #         result += r.monomial(0, i).to_free_algebra_element() @ r.monomial(0, trm_code[1] - i).to_free_algebra_element()
    #     #     return result
        
    #     # for i in range(trm_code[0] + 1):
    #     #     for j in range(trm_code[1] + 1):
    #     #         result += r.monomial(i, j).to_free_algebra_element() @r.monomial(trm_code[0] - i, trm_code[1] - j).to_free_algebra_element()
    #     # if trm_code[0] <= trm_code[1]:
    #     #     for j in range(1, trm_code[0] + 1):
    #     #         result -= coproduct(uncode([trm_code[0] -j , trm_code[1] + j]), 2)
    #     #     return result
    #     # # result -= coproduct(uncode([trm_code[0] - 1, trm_code[1] + 1]), 2)
    #     # # for i in range(trm_code[1] + 1):
    #     # #     for j in range(trm_code[0] + 1):
    #     # #         result -= r.monomial(i, j).to_free_algebra_element() @r.monomial(trm_code[1] - i, trm_code[0] - j).to_free_algebra_element()
    #     # #result -= coproduct(uncode([trm_code[1], trm_code[0]]), 2)
    #     # # for i in range(trm_code[1]):
    #     # #     for j in range(trm_code[0] + 2):
    #     # #         result -= r.monomial(i, j).to_free_algebra_element() @r.monomial(trm_code[1] - 1 - i, trm_code[0] + 1 - j).to_free_algebra_element()
    #     # for j in range(trm_code[1] + 1):
    #     #     result -= coproduct(uncode([trm_code[1] -j , trm_code[0] + j]), 2)
    #     return result
    #     # for j in range(trm_code[1] + 1):
    #     #     result -= coproduct(uncode([trm_code[0] - j, trm_code[1] + j]), 2)
    #     # return result
    #     #     result -= coproduct(uncode([trm_code[0] - 1, trm_code[1] + 1]), 2)
    #     #     for r1 in range(trm_code[0] - 1):
    #     #         for s1 in range(trm_code[1] + 3):
    #     #             result += r.monomial(r1, s1).to_free_algebra_element() @ r.monomial(trm_code[0] - 2 - r1, trm_code[1] + 2 - s1).to_free_algebra_element()
    #     #     return result
    #     # for a in range(trm_code[0] + 1):
    #     #     for b in range(trm_code[1] + 1):
    #     #         result += r.monomial(a, b).to_free_algebra_element() @ r.monomial(trm_code[0] - a, trm_code[1] - b).to_free_algebra_element()
    #     # for i in range(trm_code[1] + 1):
    #     #     for j in range(trm_code[0] + 1):
    #     #         result -= r.monomial(i, j).to_free_algebra_element() @ r.monomial(trm_code[1] - i, trm_code[0] - j).to_free_algebra_element()
    #     # result += coproduct(uncode([trm_code[1] - 1, trm_code[0] + 1]), 2)
    #     # for t in range(trm_code[1] - 1):
    #     #     for u in range(trm_code[0] + 3):
    #     #         result -= r.monomial(t, u).to_free_algebra_element() @ r.monomial(trm_code[1] - 2 - t, trm_code[0] + 2 - u).to_free_algebra_element()
    #     return result
            
    
    first_word = r.monomial(*pad_tuple(perm.trimcode, length))
    first_coprod = FA(*pad_tuple(perm.trimcode, length)).coproduct()

    first_rc_coprod = 0
    the_prin = RCGraph.principal_rc(perm, length)
    if the_weight is None:
        the_weight = the_prin.to_highest_weight()[0].length_vector

    for (word1, word2), coeff in first_coprod.items():
        the_tensors = r.monomial(*word1) @ r.monomial(*word2)
        for (rc1, rc2), coeff in the_tensors.items():
            up_prin = the_prin
            tensor = CrystalGraphTensor(rc1, rc2)
            good = True
            done = False
            while good and not done:
                done = True
                for i in range(1,len(up_prin) - 1):
                    test_raise = up_prin.raising_operator(i)
                    if test_raise is not None:
                        done = False
                        up_prin = test_raise
                        test_tensor = tensor.raising_operator(i)
                        if test_tensor is not None:
                            tensor = test_tensor
                        else:
                            good = False
                            break
                        break
            if good:
                rc11, rc22 = tensor.factors
                first_rc_coprod += coeff * r(rc11) @ r(rc22)
                #assert tensor.is_highest_weight
    #first_rc_coprod += coeff * r.from_free_algebra_element(FA(*word1)).resize(length).to_free_algebra_element() @ r.from_free_algebra_element(FA(*word2)).resize(length).to_free_algebra_element()
    if len(first_word) == 1:
        return first_rc_coprod
    for rc, coeff0 in first_word.items():
        if rc.is_principal:
            continue
        first_rc_coprod -= coeff0 * coproduct(rc.perm, length)
    #end_keys = set(first_rc_coprod.keys())
    # result_dict = {}
    # for (rc1, rc2), coeff in first_rc_coprod.items():
    #     tensor = CrystalGraphTensor(rc1, rc2)
    #     tensor_hw, _ = tensor.to_highest_weight()
    #     #if tensor_hw.crystal_weight == the_weight:
    #     rc11, rc22 = tensor_hw.factors
    #     #if tensor_hw.crystal_weight == the_weight:
    #     result_dict[(rc1.perm, rc2.perm)] = result_dict.get((rc1.perm, rc2.perm), 0) + coeff * r(rc11) @ r(rc22)

    # #assert all(v>=0 or CrystalGraphTensor(rc1,rc2).crystal_weight != the_weight for (rc1, rc2), v in first_rc_coprod.items()), f"Negative coefficient in coproduct of {perm.trimcode}, got {first_rc_coprod}"
    return first_rc_coprod
    #return sum([v for _, v in result_dict.items()])
 
if __name__ == "__main__":
    import sys
    from schubmult.rings.free_algebra import *
    CB = FreeAlgebra(CompositionSchubertBasis)
    for perm in Permutation.all_permutations(int(sys.argv[1])):
        if perm.inv == 0:
            continue
        # if len(perm.trimcode) > 2:
        #     continue
        print(f"Testing {perm.trimcode}")
        coprod = sum([coeff * (CB(*pad_tuple(rc1.perm.trimcode, len(rc1))) @ CB(*pad_tuple(rc2.perm.trimcode, len(rc2)))) for (rc1, rc2), coeff in crystal_coprod(perm, len(perm.trimcode)).items()])
        expected = CB(*perm.trimcode).coproduct()
        assert coprod.almosteq(expected), f"Failed for {perm.trimcode}, got \n{coprod}, expected \n{expected}\ndiff:\n{coprod - expected}"
        print("bench potato")