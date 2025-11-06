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
def decompose_tensor_product(dom, u_hw_rc, length, n):
    # global hw_rc_sets
    u = u_hw_rc.perm
    crystals = {}
    highest_weights = set()
    perm_set = set((Sx(u)*Sx(dom.perm)).keys())
    #for u_hw_rc in RCGraph.all_rc_graphs(u, n - 1):
    
    for w in perm_set:
        # if not u.bruhat_leq(w):
        #     continue
        # if not dom.perm.bruhat_leq(w):
        #     continue

        # print(f"Moving on to {u=} {w=} {dom.perm=}")
        # if w not in hw_rc_sets:
        #     hw_rc_sets[w] = set()
        #     for rc_w in RCGraph.all_rc_graphs(w, n - 1):
        #         # pretty_print(rc_w)
        #         hw_rc_sets[w].add(rc_w.to_highest_weight(length=length)[0])
        for rc_w in RCGraph.all_hw_rcs(w, n - 1):
            # pretty_print(rc_w)
            got_one = False
            #rc_lw = rc_w.to_lowest_weight()[0]
            high_weight = rc_w.length_vector
            for rc_lw in rc_w.full_crystal: 
                reduced_word = rc_lw.reduced_word
                # for subword in all_reduced_subwords(reduced_word, u):
                #     compatible_seq = [MarkedInteger(a) if index in subword else a for index, a in enumerate(rc_lw.compatible_sequence)]
                #     u_tab = RootTableau.root_insert_rsk(reduced_word, compatible_seq)
                #     last_inv = 1000
                #     while u_tab.perm.inv < last_inv:
                #         last_inv = u_tab.perm.inv
                #         for box in u_tab.iter_boxes:
                #             if not isinstance(u_tab[box][1], MarkedInteger):
                #                 u_tab_test = u_tab.delete_box(box)
                #                 if u_tab_test is not None:
                #                     u_tab = u_tab_test
                #                     break
                #     if u_tab.perm.inv > u.inv:
                #         # didn't make it
                #         print(f"not sure this should happen {u_tab.perm=}!!!!!!!!!!!!111")
                #         pretty_print(u_tab)
                #         input()
                #         continue

                #     u_tab = u_tab.rectify()
                #     u_hw_rc = u_tab.rc_graph.resize(n - 1)
                #     assert u_hw_rc.perm == u

                hw_checked = set()
                for u_tab2 in u_hw_rc.full_crystal:
                    tensor = CrystalGraphTensor(dom.rc_graph.resize(n - 1), u_tab2.resize(n - 1))
                    # print(f"{tensor=}")
                    tc_elem = tensor.to_highest_weight(length=length)[0]
                    # pretty_print(tc_elem)
                    if tc_elem in hw_checked:
                        # print("Already checked")
                        # print(f"{highest_weights=}")
                        continue
                    # needed!!!
                    if tc_elem in highest_weights:
                        # print("Already known highest weight mapped to some demazure crystal")
                        continue
                    u_tab_hw = tc_elem.factors[1]
                    # hw_checked.add(tc_elem)
                    #pretty_print(dom.rc_graph)
                    # CHEK TO HIGHEST WEIGHT
                    try:
                        assert tc_elem.crystal_weight == tuple([a + b for a,b in zip_longest(dom.rc_graph.length_vector, u_tab_hw.length_vector, fillvalue=0)]), f"{tc_elem.crystal_weight=} vs {tuple([a + b for a,b in zip_longest(dom.rc_graph.length_vector, u_tab_hw.length_vector, fillvalue=0)])}"
                    except AssertionError as e:
                        print(e)
                        pretty_print(tc_elem)
                        #pretty_print(u_tab2)
                        print("domrc")
                        pretty_print(dom.rc_graph)
                        print("utab")
                        pretty_print(u_tab_hw)
                        raise
                    high_weight_check = tuple([a for a, b in zip_longest(high_weight, tc_elem.crystal_weight, fillvalue=0)])
                    low_weight_check = tuple([a for a, b in zip_longest(rc_w.to_lowest_weight(length=length)[0].length_vector, tc_elem.crystal_weight, fillvalue=0)])
                    if tc_elem.crystal_weight == high_weight_check and tc_elem.to_lowest_weight(length=length)[0].crystal_weight == low_weight_check:
                        crystals[rc_w] = crystals.get(rc_w, set())
                        crystals[rc_w].add(tc_elem)
                        highest_weights.add(tc_elem)
        #                 got_one = True
        # try:
        #     assert got_one
        # except AssertionError:
        #     print("Failed to find decomposition for")
        #     print(f"{dom.perm=} and {u=}")
        #     print("Trying to find rc_w:")
        #     pretty_print(rc_w)
        #     print("with compatible sequence:")
        #     pretty_print(rc_w.compatible_sequence)
        #     raise
    return crystals

if __name__ == "__main__":
    
    ASx = FreeAlgebra(SchubertBasis)
    n = int(sys.argv[1])
    k = int(sys.argv[2])
    assert n - 1 == k

    perms = Permutation.all_permutations(n)

    hw_tabs = set()
    for perm in perms:
        if perm.minimal_dominant_above() != perm:
            continue
        hw_tabs.update([RootTableau.from_rc_graph(rc.to_highest_weight(length=None)[0]) for rc in RCGraph.all_rc_graphs(perm, n - 1)])



    rc_ring = RCGraphRing()
    tot_suc = 0
    # u times v
    print("NOTE THIS IS CORRECT AND AN ASSOCIATIVE ACTION FOR DOMINANT PERMS")
    the_schubs = {}
    # top_k_dom = {}
    used = {}
    for hw_tab0, v in itertools.product(hw_tabs, perms):
        # TEMP DOM TEST
        if hw_tab0.perm.inv == 0 or v.inv == 0:# or len(v.descents()) > 1 or k - 1 not in v.descents() or not set(v.trimcode).issubset({0,1}):
            continue
        # DOM TEST
        if hw_tab0.perm.minimal_dominant_above() != hw_tab0.perm:
            print("TEMP DOM TEST")
            continue
        hw_tab = hw_tab0
        print("hw_tab")
   
        pretty_print(hw_tab.rc_graph)
        exclude = set()
        for v_rc in RCGraph.all_rc_graphs(v, n - 1):
            # if not v_rc.is_principal:
            #     continue
            if v_rc in exclude:
                continue
            crystals = decompose_tensor_product(hw_tab, v_rc, length=None, n=n)
            
            print("Product:")
            pretty_print(hw_tab0.rc_graph)
            print("and")
            pretty_print(v_rc)
        # MUST MODIFY SM. RULE: WEIGHT PRESERVING, DIVDIFF from div_perm
        # THIS IS CRYSTAL LEVEL
        
            print("Product:")
            pretty_print(hw_tab0.rc_graph)
            print("and")
            print(f"{v.trimcode}")
            # MUST MODIFY SM. RULE: WEIGHT PRESERVING, DIVDIFF from div_perm
            # THIS IS CRYSTAL LEVEL
            sm = the_schubs.get((hw_tab0.perm, v), rc_ring.zero)
            for the_rc, st in crystals.items():
                for td in st:
                    if td.factors[1] in exclude:
                        continue
                    exclude.add(td.factors[1])
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
        for rc, coeff in sm.items():
            hw = rc.to_highest_weight()[0]
            if hw in dd:
                continue
            dd.add(hw)
            for rc0 in hw.full_crystal:
                the_sum += coeff * rc0.polyvalue(Sx.genset)
            #the_sum += coeff * Sx(rc.perm)
        prod2 = Sx(the_sum)
        print(f"{prod=}")
        print(f"{prod2=}")
        assert prod2 - prod == Sx.zero
        print(prod)
        

        