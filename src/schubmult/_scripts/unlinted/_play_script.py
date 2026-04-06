from multiprocessing import Pool, cpu_count
from schubmult.rings.polynomial_algebra import *
from schubmult import *
from schubmult.abc import x
from schubmult.utils.perm_utils import mu_A
from schubmult.rings.free_algebra import *
from schubmult.utils.tuple_utils import pad_tuple
from schubmult.combinatorics.permutation import bad_classical_patterns
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

neg_pats = [Permutation([1, 4, 3, 2]), Permutation([4, 1, 3, 2]), Permutation([3, 1, 4, 2])]
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

ring = GrassTensorAlgebra()
g = GrassTensorAlgebra()
r = RCGraphRing()
def grass_tensor_elem(perm, n = None):
    if n is None:
        n = len(perm.trimcode)
    grass_elem = ring.zero
    cem = RCGraph.full_CEM(perm, n)
    
    for rc, cem_dict in cem.items():
        #if len(rc.perm) <= n:
        if True:
            grass_elem += ring.from_dict(cem_dict)
        elif rc.perm != perm:
            grass_elem += ring.from_dict(cem_dict)
        else:
            rc_base, rc_grass = rc.resize(n).squash_decomp()
            if rc_grass.perm.inv != 0:
                grass_elem += ring.from_dict({(*k, rc_grass): v for k, v in RCGraph.full_CEM(rc_base.perm, n - 1)[rc_base.resize(n - 1)].items()})
            else:
                grass_elem += ring.from_dict(cem_dict)
    return grass_elem

def doit(grass_elem1, n):
    picnic1 = (g@r).zero
    for key, coeff in grass_elem1.items():
        to_add1 = g.one
        index = 0
        for rc in key:
            if len(rc.perm.trimcode) < n - 1:
                to_add1 *= g((rc,))
                index += 1
            else:
                break
        assert (next(iter([len(rc.perm) for rc in to_add1.to_rc_graph_ring_element()]))) <= n - 1
        picnic1 += coeff * to_add1 @ g(key[index:]).to_rc_graph_ring_element()
    return picnic1


if __name__ == "__main__":
    from schubmult.combinatorics.crystal_graph import CrystalGraphTensor
    from schubmult.abc import *
    from schubmult.rings.polynomial_algebra import *
    from sympy import pretty_print
    import sys
    n = int(sys.argv[1])
    extra = 2
    perms = [perm for perm in Permutation.all_permutations(n + extra) if len(perm.trimcode) <= n]
    
    

    grass_tensor_elems = {}
    hw_grass = {}

    for perm1, perm2 in itertools.product(perms, repeat=2):
        if perm1.inv == 0 and perm2.inv == 0:
            continue
        if perm1.inv == 0 or perm2.inv == 0:
            continue
        cem1 = RCGraph.full_CEM(perm1, n - 1)
        elem1 = g.zero
        for rc, cem_dict in cem1.items():
            for rc_tup, coeff in cem_dict.items():
                if CrystalGraphTensor(*rc_tup).is_highest_weight:
                    elem1 += coeff * g(rc_tup)
        #pretty_print(elem1)
        elem2 = g.zero
        cem2 = RCGraph.full_CEM(perm2, n - 1)
        for rc, cem_dict in cem2.items():
            for rc_tup, coeff in cem_dict.items():
                if CrystalGraphTensor(*rc_tup).is_highest_weight:
                    elem2 += coeff * g(rc_tup)
        # assert all(v > 0 for v in elem1.values()), f"Unexpected negative coefficient for {perm1=}"
        # assert all(v > 0 for v in elem2.values()), f"Unexpected negative coefficient for {perm2=}"
        result = (r@r).zero
        prd = Sx(perm1) * Sx(perm2)
        for key1, coeff1 in elem1.items():
            for key2, coeff2 in elem2.items():
                for key11 in key1.full_crystal:
                    for key22 in key2.full_crystal:
                        tensor1 = next(iter(doit(g(key11), n + 1)))
                        tensor2 = next(iter(doit(g(key22), n + 1)))
                        other_tensor = next(iter(g(tensor1[0]) * g(tensor2[0])))
                        if other_tensor.is_highest_weight:
                            base_tensor = next(iter(doit(g(other_tensor), n)))
                            second_rc = next(iter((g((base_tensor[1],)) * g((tensor1[1],)) * g((tensor2[1],))).to_rc_graph_ring_element()))
                            tensor2 = CrystalGraphTensor(*base_tensor[0].factors, second_rc)
                            if tensor2.is_highest_weight:
                                result += coeff1 * coeff2 * g(base_tensor[0]).to_rc_graph_ring_element() @ r(second_rc)

                        
                        # if CrystalGraphTensor(next(iter(g(tensor[0]).to_rc_graph_ring_element())),tensor[1]).is_highest_weight:
                        #     result += coeff1 * coeff2 * g(tensor[0]).to_rc_graph_ring_element() @ r(tensor[1])
                        # tensor2 = next(iter(doit(g(key22), n)))
                        # if len(tensor1[0].factors) == 0 and len(tensor2[0].factors) == 0:
                        #     result += coeff1 * coeff2 * r(RCGraph([])) @ r(tensor1[1].resize(len(tensor2[1])).squash_product(tensor2[1]))
                        # else:
                        #     if len(tensor1[0].factors) == 0:
                        #         tensor = tensor2[0]
                        #     elif len(tensor2[0].factors) == 0:
                        #         tensor = tensor1[0]
                        #     else:
                        #         tensor = CrystalGraphTensor(*tensor1[0].factors, *tensor2[0].factors)
                        #     if tensor.is_highest_weight:
                        #         #assert coeff1 * coeff2 > 0, f"Unexpected negative coefficient for {perm1=}, {perm2=}"
                        #         # rc_elem = g(tensor).to_rc_graph_ring_element()
                        #         # #tensorbob = CrystalGraphTensor(rc_elem, tensor1[1].resize(len(tensor2[1])).squash_product(tensor2[1]))
                                
                        #         # if tensor2[1].perm.inv !=0 and tensor1[1].perm.inv != 0:
                        #         #     tensorbob = CrystalGraphTensor(rc_elem, r(tensor1[1].resize(len(tensor2[1])).squash_product(tensor2[1])))
                        #         # elif tensor2[1].perm.inv != 0:
                        #         #     tensorbob = CrystalGraphTensor(rc_elem, r(tensor2[1]))
                        #         # elif tensor1[1].perm.inv != 0:
                        #         #     tensorbob = CrystalGraphTensor(rc_elem, r(tensor1[1]))
                        #         # else:
                        #         #     tensorbob = CrystalGraphTensor(rc_elem, r(RCGraph([])).resize(len(rc_elem)))
                        #         # if tensorbob.is_highest_weight:
                        #         #     result += coeff1 * coeff2 * (r@4)(tensorbob.factors)
                        #         rc1 = next(iter(g(tensor).to_rc_graph_ring_element()))
                        #         rc2 = next(iter((g((tensor1[1],)) * g((tensor2[1],))).to_rc_graph_ring_element()))
                        #         tensor2 = CrystalGraphTensor(rc1, rc2)
                        #         if tensor2.is_highest_weight:
                        #             result += coeff1 * coeff2 * (r@r)(tensor2.factors)



                                    
                # tensor = next(iter(g(key1) * g(key2)))
                # if tensor.is_highest_weight:
                #     result += coeff1 * coeff2 * g(tensor).to_rc_graph_ring_element()
        
        prd2 = Sx.zero
        for (rc_base, rc_grass), coeff in result.items():
            if rc_base.perm.inv == 0:
                rc = rc_grass
            elif rc_grass.perm.inv == 0:
                rc = rc_base
            else:
                rc = rc_base.resize(len(rc_grass)).squash_product(rc_grass)
            assert prd.get(rc.perm, 0) ==coeff, f"Failure for {perm1}, {perm2}, got {rc.perm}, {prd=}, {rc=}, {coeff=}"
            if rc.perm not in prd2:
                prd2 += coeff * Sx(rc.perm)
        assert prd == prd2, f"Failure for {perm1}, {perm2}, got {prd=}, {prd2=}"
        print("Success for", perm1, perm2)


    # for perm in perms:
    #     if perm.inv == 0:
    #         continue
    #     cem = RCGraph.full_CEM(perm, n - 1)
    #     elem = r.zero
    #     for rc, cem_dict in cem.items():
    #         for rc_tup, coeff in cem_dict.items():
    #             if CrystalGraphTensor(*rc_tup).is_highest_weight:
    #                 elem += coeff * g(rc_tup).to_rc_graph_ring_element()
    #     pretty_print(elem)

    # for perm1, perm2 in itertools.product(perms, repeat=2):
    #     if perm1.inv == 0 or perm2.inv == 0:
    #         continue
    #     cems1 = RCGraph.full_CEM(perm1, n - 1)
    #     hw_set = set()
    #     cem_result1 = g.zero
    #     for rc, cem_dict in cems1.items():
    #         # if rc.perm != perm1:
    #         #     continue
    #         cem_result = g.zero
    #         for rc_tup, coeff in cem_dict.items():
    #             cem_result1 += coeff * g(rc_tup)
        
    #     cems2 = RCGraph.full_CEM(perm2, n - 1)
    #     hw_set = set()
    #     cem_result2 = g.zero
    #     for rc, cem_dict in cems2.items():
    #         # if rc.perm != perm2:
    #         #     continue
    #         for rc_tup, coeff in cem_dict.items():
    #             cem_result2 += coeff * g(rc_tup)
    #     result = r.zero
    #     for key1, coeff1 in cem_result1.items():
    #         for key2, coeff2 in cem_result2.items():
    #             # first_result = next(iter((g(key1[0]) * g(key2[0])).to_rc_graph_ring_element())).resize(n)
    #             # tensor = CrystalGraphTensor(first_result.resize(n),key1[1].resize(n).squash_product(key2[1].resize(n))).to_highest_weight()[0]
    #             # result += coeff1 * coeff2 * r(tensor.factors[0].squash_product(tensor.factors[1]).reverse_raise_seq(tensor.raise_seq))
    #             # assert coeff1 > 0
    #             # assert coeff2 > 0
    #             # rc0 = next(iter((g(key1))))
    #             # rc1 = next(iter(g(key2)))#next(iter(g(key2).to_rc_graph_ring_element()))
    #             key3 = next(iter(g(key1) *g(key2)))
    #             if key3.is_highest_weight:
    #                 rc_to_add = g(key3).to_rc_graph_ring_element()
    #                 result += coeff1 * coeff2 * rc_to_add
    #             # for rc01 in rc0.full_crystal:
    #             #     for rc11 in rc1.full_crystal:
    #             #         spacious = next(iter(g(rc01) * g(rc11)))
    #             #         if spacious.is_highest_weight:
    #             #             rc_to_add = g(spacious).to_rc_graph_ring_element()
    #             #             result += coeff1 * coeff2 * rc_to_add
    #     prd = Sx(perm1) * Sx(perm2)
    #     prd2 = Sx.zero
    #     for rc, coeff in result.items():
    #         assert prd.get(rc.perm, 0) ==coeff, f"Failure for {perm1}, {perm2}, got {rc.perm}, {prd=}, {rc=}"
    #         if rc.perm not in prd2:
    #             prd2 += coeff * Sx(rc.perm)
    #     assert prd == prd2, f"Failure for {perm1}, {perm2}, got {prd=}, {prd2=}"
    #     print("Success for", perm1, perm2)
    # #     # rcs = [rc for rc in cems if rc.perm != perm]
    # #     # for rc in rcs:
    # #     #     test_cem = cems[rc]
    # #     #     assert ring.from_dict(test_cem).to_rc_graph_ring_element().almosteq(r.zero)
    # #     #     for k in range(1, n):
    # #     #         for p in range(1, k + 1):
    # #     #             for try_rc in RCGraph.all_rc_graphs(uncode([1]*k), k):
    # #     #                 sumup = ring.zero
    # #     #                 for cemm, coeff in test_cem.items():
    # #     #                     sumup += coeff * ring(cemm) * ring((try_rc,))
    # #     #                 if sumup.to_rc_graph_ring_element() != r.zero:
    # #     #                     print(f"Failed on {perm} with k={k}, got {sumup}")
    # #     #                     print(f"Tested against {try_rc}")
    # #     #                     print(f"CEM was {test_cem}")
    # #     #                     print(f"RC was {rc}")
    # #     #                     input()
    # #     #                 left_sumup = ring.zero
    # #     #                 for cemm, coeff in test_cem.items():
    # #     #                     left_sumup += coeff * ring((try_rc,)) * ring(cemm)
    # #     #                 if left_sumup.to_rc_graph_ring_element() != r.zero:
    # #     #                     print(f"Failed on {perm} with k={k}, got {left_sumup}")
    # #     #                     print(f"Tested against {try_rc}")
    # #     #                     print(f"CEM was {test_cem}")
    # #     #                     print(f"RC was {rc}")
    # #     #                     input()

    # #     grass_tensor_elems[perm] = ring.zero
    # #     hw_grass[perm] = ring.zero
    # #     #part = tuple(Permutation.w0(n).trimcode)
    # #     part = None
    # #     cem = RCGraph.full_CEM(perm, n - 1, part)
    # #     for rc, cem_dict in cem.items():
    # #         valval = ring.from_dict(cem_dict)
    # #         #if True:
    # #         # if len(rc.perm) > n:
    # #         #     continue
    # #         grass_tensor_elems[perm] += valval
    # #     # if perm.inv > 0:
    # #     #     hw_keys = set()
    # #     #     for rc in grass_tensor_elems[perm].keys():
    # #     #         if len(rc) == 1:
    # #     #             hw_keys.add((rc[0].to_highest_weight()[0],))
    # #     #         else:
    # #     #             hw_keys.add(CrystalGraphTensor(*rc).to_highest_weight()[0].factors)
    # #         # hw_grass[perm] = sum([coeff * ring(key) for key, coeff in grass_tensor_elems[perm].items() if key in hw_keys])
    # #         # if any(coeff < 0 for key, coeff in hw_grass[perm].items()):
    # #         #     print(f"Neg coeff {perm=}")
    # #         pretty_print(grass_tensor_elems[perm])
    # #         #     assert any(perm.has_pattern(pat) for pat in neg_pats), f"Unexpected negative coefficient for {perm=}"
    # #         # else:
    # #         #     assert all(not perm.has_pattern(pat) for pat in neg_pats), f"Unexpected negative coefficient for {perm=}"
                
    #input()
    sys.exit(0)
    extra = 3
    codeperms = [perm for perm in Permutation.all_permutations(n + extra) if len(perm.trimcode) < n]
    for perm1, perm2 in itertools.product(codeperms, repeat=2):
        # if perm1 != uncode([0,0,1]) or perm2 != Permutation([1,4,3,2]):
        #     continue
        # if any(perm1.has_pattern(pat) for pat in neg_pats):# or any(perm2.has_pattern(pat) for pat in neg_pats):
        #     continue
        # prd = Sx(perm1) * Sx(perm2)
        # if len(perm1) <= n and len(perm2) <= n:
        #     continue
        grass_elem1 = grass_tensor_elem(perm1, n - 1)
        grass_elem2 = grass_tensor_elem(perm2, n - 1)
        picnic1 = (ring @ r).zero
        picnic2 = (ring @ r).zero
        for key, coeff in grass_elem1.items():
            to_add1 = ring.one
            index = 0
            for rc in key:
                if len(rc.perm.trimcode) < n - 1:
                    to_add1 *= ring((rc,))
                    index += 1
                else:
                    break
            assert (next(iter([len(rc.perm) for rc in to_add1.to_rc_graph_ring_element()]))) <= n
            picnic1 += coeff * to_add1 @ ring(key[index:]).to_rc_graph_ring_element()
        # pretty_print(picnic1)
        #assert all(v > 0 for v in picnic1.values()), f"Unexpected negative coefficient for {perm1=}"
        for key, coeff in grass_elem2.items():
            to_add2 = ring.one
            index = 0
            for rc in key:
                if len(rc.perm.trimcode) < n - 1:
                    to_add2 *= ring((rc,))
                    index += 1
                else:
                    break
            assert (next(iter([len(rc.perm) for rc in to_add2.to_rc_graph_ring_element()]))) <= n
            picnic2 += coeff * to_add2 @ ring(key[index:]).to_rc_graph_ring_element()
        potato = r.zero
        for key1, coeff1 in picnic1.items():
            for key2, coeff2 in picnic2.items():
                potato0 = coeff1 * coeff2 * (ring(key1[0]) * ring(key2[0])).to_rc_graph_ring_element()
                rc1 = key1[1].resize(n - 1)
                rc2 = key2[1].resize(n - 1)
                for rcc, coeff in potato0.items():
                    potato += coeff * r(rcc.resize(n - 1).squash_product(rc1).squash_product(rc2))
        prd_rc = potato
        prd = Sx(perm1) * Sx(perm2)
        #pretty_print(picnic2)
        for perm0, coeff in prd.items():
            # if rc.is_highest_weight and rc.resize(n).squash_decomp()[1].perm.inv == 0:
            # #if len(rc.perm) <= n:
            assert prd_rc.get(RCGraph.principal_rc(perm0, n - 1), 0) == coeff, f"Failed on {perm1} * {perm2}, got {rc.perm}: {coeff} which is not in {prd}\n{prd.get(rc.perm,0)=}\n{rc.perm=}"
        print(f"Checked {perm1} * {perm2}")
    #input()
    # for perm1, perm2 in itertools.product(codeperms, repeat=2):
    #     # if perm1 != uncode([0,0,1]) or perm2 != Permutation([1,4,3,2]):
    #     #     continue
    #     # if any(perm1.has_pattern(pat) for pat in neg_pats):# or any(perm2.has_pattern(pat) for pat in neg_pats):
    #     #     continue
    #     prd = Sx(perm1) * Sx(perm2)
    #     # if len(perm1) <= n and len(perm2) <= n:
    #     #     continue
    #     prd_elem = ring.zero
    #     grass_elem1 = grass_tensor_elem(perm1, n)
    #     grass_elem2 = grass_tensor_elem(perm2, n)
    #     # for rc1, cem_dct in RCGraph.full_CEM(perm1, n -1):
    #     #     rc1_split, rc1_grass = rc1.resize(n).squash_decomp()
    #     #     grass_elem_base = RCGraph.full_CEM(rc1_split.perm, n -1)[rc1_split.resize(n - 1)]
    #     #     new_dct = {(k, rc1_grass): v for k, v in grass_elem_base.items()}
    #     #     grass_elem1 += ring.from_dict(new_dct)

    #     # for rc2, coeff2 in RCGraph.full_CEM(perm2, n - 1):
    #     #     rc2_split, rc2_grass = rc2.resize(n).squash_decomp()
                
    #     prd_elem = grass_elem1 * grass_elem2
    #     #pretty_print(prd_elem)

    #     prd_rc = prd_elem.to_rc_graph_ring_element()
    #     #atleastone = False
    #     for rc, coeff in prd_rc.items():
    #         # if rc.is_highest_weight and rc.resize(n).squash_decomp()[1].perm.inv == 0:
    #         # #if len(rc.perm) <= n:
    #         assert prd.get(rc.perm, 0) == coeff, f"Failed on {perm1} * {perm2}, got {rc.perm}: {coeff} which is not in {prd}\n{prd.get(rc.perm,0)=}\n{rc.perm=}"
    #             #atleastone = True
                
        #assert atleastone or all(len(p) > n for p in prd), f"Failed on {perm1} * {perm2}, got {prd_rc} which has no highest weight terms with trivial right factor, expected {prd}"
        #print(f"Checked {perm1} * {perm2}")