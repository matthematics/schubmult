import copy
import itertools
import sys
from functools import cache
from itertools import zip_longest

import sympy
from sympy import pretty_print

from schubmult import CrystalGraphTensor, FreeAlgebra, Permutation, RCGraph, RCGraphRing, RootTableau, SchubertBasis, Sx


def all_reduced_subwords(reduced_word, u):
    if u.inv > len(reduced_word):
        return set()
    if u.inv == 0:
        return {()}
    ret_set = set()
    for index in range(len(reduced_word) - 1, -1, -1):
        a = reduced_word[index]
        if a - 1 in u.descents():
            new_u = u.swap(a - 1, a)
            old_set = all_reduced_subwords(reduced_word[:index], new_u)
            for subword in old_set:
                new_subword = (*subword, index)
                ret_set.add(new_subword)
    return ret_set

class MarkedInteger(int):
    pass

hw_rc_sets = {}
def decompose_tensor_product(dom, u_rc, n):
    from schubmult import RCGraphRing
    rc_ring = RCGraphRing()
    tring = rc_ring @ rc_ring
    # global hw_rc_sets
    crystals = {}
    assert len(u_rc) == n - 1
    assert len(dom) == n - 1
    if u_rc.inv == 0:
        crystals[dom] = {CrystalGraphTensor(dom, u_rc)}
        return crystals
    if dom.inv == 0:
        crystals[u_rc] = {CrystalGraphTensor(dom, u_rc)}
        return crystals
    if len(u_rc) == 0:
        assert len(dom) == 0
        crystals[dom] = {CrystalGraphTensor(dom, u_rc)}
        return crystals
    if len(u_rc) == 1:
        assert len(dom) == 1
        crystals[RCGraph.one_row(len(dom[0]) + len(u_rc[0]))] = {CrystalGraphTensor(dom, u_rc)}
        return crystals
    cut_dom = dom.vertical_cut(n-2)[0]
    cut_u = u_rc.vertical_cut(n-2)[0]
    # print("Cutting:")
    # pretty_print(cut_dom)
    # pretty_print(cut_u)
    cut_crystals = decompose_tensor_product(cut_dom, cut_u, n - 1)
    # print(f"{cut_crystals=}")

    for rc_w_cut, tensor_elems in cut_crystals.items():
        up_rc =  rc_ring(rc_w_cut) * rc_ring(RCGraph.one_row(len(dom[-1]) + len(u_rc[-1])))
        up_tensor = tring.zero
        for t_elem in tensor_elems:
            to_add =  tring(t_elem.factors) * tring((RCGraph.one_row(len(dom[-1])),RCGraph.one_row(len(u_rc[-1]))))
            # pretty_print(to_add)
            for (rc1, rc2), coeff in to_add.items():
                assert coeff == 1
                if rc1 != dom or rc2 != u_rc:
                    continue
                # tcryst = CrystalGraphTensor(rc1, rc2)
                # for tw in tcryst.full_crystal:
                #     if tw.factors[1] == u_rc:
                up_tensor += coeff * tring((rc1, rc2))
                #        break
                    #up_tensor += coeff * tring((rc1, rc2))
        # print("up_tensor=")
        # pretty_print(up_tensor)
        # print("up_rc=")
        # pretty_print(up_rc)
        for w_rc, coeff in up_rc.items():
            assert coeff == 1
            high_weight = w_rc.to_highest_weight()[0].crystal_weight
            low_weight = w_rc.to_lowest_weight()[0].crystal_weight
            if w_rc.perm not in (Sx(dom.perm) * Sx(u_rc.perm)).keys():
                continue
            for (rc1, u_rc2), coeff2 in up_tensor.items():
                assert coeff2 == 1
                tensor = CrystalGraphTensor(rc1, u_rc2)
                tensor_hw = tensor.to_highest_weight()[0]
                tensor_lw = tensor.to_lowest_weight()[0]
                if tensor_hw.crystal_weight == high_weight and tensor_lw.crystal_weight == low_weight:
                    crystals[w_rc] = crystals.get(w_rc, set())
                    crystals[w_rc].add(tensor)
    try:
        assert len(crystals) == 1
    except AssertionError:
        print("Error: More than one crystal found.")
        for rc, tensor_set in crystals.items():
            pretty_print(rc)
            for tensor in tensor_set:
                pretty_print(tensor)
                print("high weight")
                pretty_print(rc.to_highest_weight()[0])
                pretty_print(tensor.to_highest_weight()[0])
                print("low weight")
                pretty_print(rc.to_lowest_weight()[0])
                pretty_print(tensor.to_lowest_weight()[0])
        raise
    return crystals

def is_decomposable(w):
    for i in range(1, len(w) - 1):
        coset, w_J = w.coset_decomp(*list(range(1, i + 1)),*list(range(i + 2, len(w)))) 
        if coset.inv == 0 and set(w_J.code[:i+1]) != {0} and set(w_J.code[i+2:]) != {0}:
            return True
    return False


if __name__ == "__main__":
    
    ASx = FreeAlgebra(SchubertBasis)
    n = int(sys.argv[1])
    # k = int(sys.argv[2])
    # assert n - 1 == k

    perms = Permutation.all_permutations(n)

    hw_tabs = set()
    for perm in perms:
        if perm.minimal_dominant_above() != perm:# or perm != Permutation.ref_product(3,2,1,2):
            continue
        hw_tabs.update([rc.to_highest_weight(length=None)[0] for rc in RCGraph.all_rc_graphs(perm, n - 1)])



    rc_ring = RCGraphRing()
    tot_suc = 0
    # u times v
    #print("NOTE THIS IS CORRECT AND AN ASSOCIATIVE ACTION FOR DOMINANT PERMS")
    the_schubs = {}
    # top_k_dom = {}
    used = {}
    for hw_tab0, v in itertools.product(hw_tabs, perms):
        # TEMP DOM TEST
        if hw_tab0.perm.inv == 0 or v.inv == 0:# or len(v.descents()) > 1 or k - 1 not in v.descents() or not set(v.trimcode).issubset({0,1}):
            continue
        if v == Permutation.ref_product(2,4):
            print("SKIPPING DECOMPOSABLE")
            assert is_decomposable(v)
            continue
    
        if is_decomposable(v):
            print("SKIPPING DECOMPOSABLE")
            continue
        # DOM TEST
        if hw_tab0.perm.minimal_dominant_above() != hw_tab0.perm:
            print("TEMP DOM TEST")
            continue
        hw_tab = hw_tab0
        print("hw_tab")
   
        pretty_print(hw_tab)
        #exclude = set()
        for v_rc in RCGraph.all_rc_graphs(v, n - 1):
            print("Product:")
            pretty_print(hw_tab0)
            print("and")
            print(f"{v.trimcode}")
         
            crystals = decompose_tensor_product(hw_tab, v_rc, n)
            for rc_w in crystals:
                pretty_print(rc_w)
                pretty_print(crystals[rc_w])
         
            # MUST MODIFY SM. RULE: WEIGHT PRESERVING, DIVDIFF from div_perm
            # THIS IS CRYSTAL LEVEL
            sm = the_schubs.get((hw_tab0.perm, v), rc_ring.zero)
            for the_rc, st in crystals.items():
                for td in st:
                    sm += rc_ring.from_dict({the_rc: 1})

            pretty_print(sm)
            the_schubs[(hw_tab0.perm, v)] = sm
        u = hw_tab0.perm
        Permutation.print_as_code = True
    #for (u, v), val in the_schubs.items():
        print(f"{u} * {v}=")
        pretty_print(sm)
        prod = Sx(u) * Sx(v)
        the_sum = sympy.S.Zero
        dd = set()
        the_sum2 = Sx.zero
        for rc, coeff in sm.items():
            #hw = rc.to_highest_weight()[0]
            # if hw in dd:
            #     continue
            # dd.add(hw)
            #for rc0 in hw.full_crystal:
            the_sum += coeff * rc.polyvalue(Sx.genset)
            if rc.is_principal:
                the_sum2 += coeff * Sx(rc.perm)
            #the_sum += coeff * Sx(rc.perm)
        prod2 = Sx(the_sum)
        print(f"{prod=}")
        print(f"{prod2=}")
        print(f"{the_sum2=}")
        assert prod2 - prod == Sx.zero
        assert prod == the_sum2
        print(prod)
        

        