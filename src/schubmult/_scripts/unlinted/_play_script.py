from multiprocessing import Pool, cpu_count
from schubmult.rings.polynomial_algebra import *
from schubmult import *
from schubmult.abc import x
from schubmult.utils.perm_utils import mu_A
from schubmult.rings.free_algebra import *
from schubmult.utils.tuple_utils import pad_tuple
import itertools

def check_perm(perm_arr):
    perm = Permutation(perm_arr)
    n = len(perm_arr)
    r = RCGraphRing()
    for i in range(len(perm.trimcode), n):
        rc_elem = r.from_free_algebra_element(ASx(perm, i))
        if not all(int(v) in {1, 0, -1} for v in rc_elem.values()):
            return (False, perm_arr, dict(rc_elem))
    return (True, perm_arr, None)

EP = PolynomialAlgebra(ElemSymPolyBasis(Sx.genset))

# def in_good_algebra(elem_comp, monomial_vec):
#     #monom = EP.from_dict({elem_comp: 1}).change_basis(EP._basis.monomial_basis)
#     #return monom.get(monomial_vec, 0) != 0
#     assert isinstance(elem_comp[1], Permutation)
#     if FA(*)
#     return True

SS = FreeAlgebra(SchubertSchurBasis)
EE = FreeAlgebra(ElementaryBasis)

def convert_tensor(tensor):
        ret = (ASx@FA@ASx@FA).zero
        for (rc1, rc2), coeff in tensor.items():
            ret += coeff * ASx(rc1.perm, len(rc1)) @ FA(*rc1.length_vector) @ ASx(rc2.perm, len(rc2)) @ FA(*rc2.length_vector)
        return ret

def schub_weight_equal(tensor1, tensor2):
    
    
    return convert_tensor(tensor1).almosteq(convert_tensor(tensor2))


r = RCGraphRing()
def rc_to_ss(rc1):
    ret = EE.zero
    if isinstance(rc1, RCGraph):
        schuby = r(rc1).to_free_algebra_element().change_basis(ElementaryBasis)
        wordy = FA(*rc1.length_vector).change_basis(ElementaryBasis)
        
        for elem_comp, coeff in schuby.items():
            if coeff == 0:
                continue
            if wordy.get(elem_comp, 0) != 0:
                ret += coeff * EE(*elem_comp)
    else:
        for rc, coeff in rc1.items():
            ret += coeff * rc_to_ss(rc)
    return ret

def full_coprod(perm, weight, length):
    cem = RCGraph.full_CEM(perm, length)
    result = 0
    for rc, cem_dict in cem.items():
        if rc.length_vector != weight:
            continue
        for rc_tup, coeff in cem_dict.items():
            if coeff == 0:
                continue
            if tuple([a + b for a, b in zip(rc_tup[0].length_vector, rc_tup[1].length_vector)]) == weight:
                yield (rc_tup, coeff)

if __name__ == "__main__":
    from schubmult.combinatorics.crystal_graph import CrystalGraphTensor
    from schubmult.abc import *
    from schubmult.rings.polynomial_algebra import *
    from sympy import pretty_print
    n = 5
    
    perms = Permutation.all_permutations(n)
    Sy = SingleSchubertRing(y)
    #for perm in perms:
        # def or_cp(p):
        #     re = 0
        #     for d in range(p + 1):
        #         re += r.monomial(d) @ r.monomial(p - d)
        #     return re
        
        # def wr_cpd(*wrd):
        #     if len(wrd) == 1:
        #         return or_cp(wrd[0])
        #     return or_cp(wrd[0]) * wr_cpd(*wrd[1:])
        # elem = ASx(perm, n - 1).coproduct()#change_basis(WordBasis)
        # # full_pudge = 0
        # cem = RCGraph.full_CEM(perm, n -1)
        # pudge = 0
        # for rc, cem_dict in cem.items():#RCGraph.all_rc_graphs(perm, n - 1):
        # #     pudge = RCGraph.multiply_reps(cem_dict, {(): 1})
        # #     for rcc, coeff in pudge.items():
        # #         full_pudge += r(rcc) @ Sy.from_expr(coeff)
        #     for ((perm1, _), (perm2, _)), coeff in elem.items():
        #             cem1 = RCGraph.full_CEM(perm1, n - 1)
        #             cem2 = RCGraph.full_CEM(perm2, n -1)    
                    
        #             for rc1, rc2 in itertools.product(cem1.keys(), cem2.keys()):
        #                 #lv = tuple([a + b for a, b in zip(rc1.resize(n-1).length_vector, rc2.resize(n-1).length_vector)])
        #                 if RCGraph.multiply_reps(cem1[rc1], cem2[rc2]).resize(n-1).almosteq(r(rc)):
        #                     pudge += RCGraph.multiply_reps(cem1[rc1], {(): 1}).resize(n-1)@RCGraph.multiply_reps(cem2[rc2], {(): 1}).resize(n-1)
                    
        #     if pudge != 0 and pudge != (r@r).zero:
        #         pretty_print(rc)
        #         pretty_print(pudge)
        # # for rcc, coeff in pudge.items():
        # #     print(Sy.from_expr(coeff))
        # #pretty_print(full_pudge)
            

    #print(Schub((uncode([0,2,1]),3)).change_basis(ElemSymPolyBasis))
    cprd = {}
    for perm1, perm2 in itertools.product(perms, repeat=2):
        # if perm1.inv == 0 or perm2.inv == 0:
        #     continue
        result = 0
        
        prd = Sx(perm1) * Sx(perm2)
        prd_rc = r.zero
        for perm, coeff in prd.items():
            prd_rc += coeff * r.schub(perm, n - 1)
    #     #perms_seen = {}
    #     seen = set()
    #     prd = Sx(perm1) * Sx(perm2)
    #     N = max(len(perm1), len(perm2))
    #     hws = set()
        cem1 = RCGraph.full_CEM(perm1, n -1)
        cem2 = RCGraph.full_CEM(perm2, n - 1)
        result = r.zero
        for rc1, rc2 in itertools.product(cem1.keys(), cem2.keys()):
            # tensor = CrystalGraphTensor(rc1_0, rc2_0).to_lowest_weight()[0]
            # if tensor in hws:
            #     continue
            # hws.add(tensor)
            # rc1, rc2 = tensor.factors   
            rc_prod = RCGraph.multiply_reps(cem1[rc1], cem2[rc2]).resize(n-1)
            # for rc, coeff in rc_prod.items():
            #     result += coeff * ASx(rc.perm, len(rc)) @ PA(*rc.length_vector).change_basis(SchubertPolyBasis)
            result += rc_prod
            for rc, fatbat in rc_prod.items():
                cprd[rc] = cprd.get(rc, (r@r).zero) + fatbat * r(rc1)@r(rc2)
        # for ((perm1, _), (perm2, _)), coeff in result.items():
        #     assert perm1 == perm2
        #     assert prd.get(perm1, 0) == coeff, f"Failed on {perm1} * {perm2}, got {result}, expected {prd}"
        for rc, coeff in result.items():
            assert prd_rc.get(rc, 0) == coeff, f"Failed on {perm1} * {perm2}, got {rc.perm}: {coeff} which is not in {prd_rc}\n{prd_rc.get(rc,0)=}\n{rc=}"
        

        print(result)
            # #result += Schub.from_dict(r.monomial(*wt).to_free_algebra_element())
            # rc3 = convert_tensor((r(rc1) * r(rc2)).coproduct())
            # test_rc3 = convert_tensor(r(rc1).coproduct() * r(rc2).coproduct())
            # assert rc3.almosteq(test_rc3), f"Fail \n{rc1}\n{rc2}\ngot\n{test_rc3 - rc3}\nwrong"
        print(f"Success {perm1}, {perm2}")
    for rc, coeff in cprd.items():
        if len(rc.perm) > n:
            continue
        print("Coproduct of ")
        pretty_print(rc)
        print("is")
        pretty_print(coeff)
    #         # assert len(rc3) == 1
    #         # assert next(iter(rc3.values())) == 1
    #         # assert all(v > 0 for v in rc3.values())
    #         #rc3 = next(iter(rc3.keys()))
            
    #         result += rc3
    #     # try:
    #     for rc3, coeff in result.items():
    #         # if coeff == 0:
    #         #     continue
    #         assert prd.get(rc3.perm, 0) == coeff, f"Failed on {perm1} * {perm2}, got {rc3.perm}: {coeff} which is not in {prd}"
    #         print("hoor")
    #     # except AssertionError as e:
    #     #     if len(rc3.perm) > N:
    #     #         continue
    #     #     print(f"Failed on {perm1} * {perm2}, got {rc3.perm} which is not in {prd}")
    #     #     raise e
    #     print("ok")
                
    #                 #perms_seen[wt] = perms_seen.get(wt, 0)  + Sx(w)
            
    #     #assert result == prd,  f"Failed on {perm1} * {perm2}, got {result}, expected {prd}"
    #     # if perm.inv == 0:
    #     #     continue
    #     # for rc in RCGraph.all_rc_graphs(perm):
    #     #     free_schub_elem = r(rc).to_free_algebra_element().change_basis(SchubertSchurBasis)
    #     #     accum_free = 0
    #     #     for elem_comp, coeff in free_schub_elem.items():
    #     #         #print(elem_comp)s
    #     #         # elem_monoms = EP.from_dict({elem_comp: 1}).change_basis(EP._basis.monomial_basis)
    #     #         # if elem_monoms.get(rc.length_vector, 0) != 0:
    #     #         #     accum_free += coeff * EE(*elem_comp)
    #     #         if in_good_algebra(elem_comp, rc.length_vector):
    #     #             accum_free += SS(*elem_comp)

    #     #     tensor_elem = accum_free @ FA(*rc.length_vector)
    #     #     new_tensor_elem = 0
    #     #     for (elem_comp1, monom_comp1, elem_comp2, monom_comp2), coeff in tensor_elem.coproduct().items():
    #     #         if coeff != 0:
    #     #             # if in_good_algebra(elem_comp1, monom_comp1) and in_good_algebra(elem_comp2, monom_comp2):
    #     #             #     new_tensor_elem += coeff * (EE(*elem_comp1) @ FA(*monom_comp1)) * (EE(*elem_comp2) @ FA(*monom_comp2))

    #     #             if sum(elem_comp1[0]) == sum(monom_comp1):
    #     #                 assert in_good_algebra(elem_comp1, monom_comp1)

    #     #             if sum(elem_comp2[0]) == sum(monom_comp2):
    #     #                 assert in_good_algebra(elem_comp2, monom_comp2)
    #     #     print(new_tensor_elem)
