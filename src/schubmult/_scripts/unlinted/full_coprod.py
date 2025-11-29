import copy
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
@cache
def decompose_tensor_product(dom, u, n):
    # global hw_rc_sets
    crystals = {}
    highest_weights = set()
    perm_set = set((Sx(u)*Sx(dom.perm)).keys())
    for w in perm_set:
        if len(w) > n:
            continue
        # if not u.bruhat_leq(w):
        #     continue
        # if not dom.perm.bruhat_leq(w):
        #     continue

        # print(f"Moving on to {u=} {w=} {dom.perm=}")
        if w not in hw_rc_sets:
            hw_rc_sets[w] = set()
            for rc_w in RCGraph.all_rc_graphs(w, n - 1):
                # pretty_print(rc_w)
                if not rc_w.is_highest_weight:
                    continue
                hw_rc_sets[w].add(rc_w)
        for rc_w in hw_rc_sets[w]:
            # pretty_print(rc_w)
            high_weight = rc_w.length_vector
            reduced_word = rc_w.reduced_word
            for subword in all_reduced_subwords(reduced_word, u):
                compatible_seq = [MarkedInteger(a) if index in subword else a for index, a in enumerate(rc_w.compatible_sequence)]
                u_tab = RootTableau.root_insert_rsk(reduced_word, compatible_seq)
                last_inv = 1000
                while u_tab.perm.inv < last_inv:
                    last_inv = u_tab.perm.inv
                    for box in u_tab.iter_boxes:
                        if not isinstance(u_tab[box][1], MarkedInteger):
                            u_tab_test = u_tab.delete_box(box)
                            if u_tab_test is not None:
                                u_tab = u_tab_test
                                break
                if u_tab.perm.inv > u.inv:
                    # didn't make it
                    continue

                u_tab = u_tab.rectify()
                u_hw_rc = u_tab.rc_graph.resize(n - 1)
                assert u_hw_rc.perm == u

                hw_checked = set()
                for u_tab2 in u_hw_rc.full_crystal:
                    tensor = CrystalGraphTensor(dom.rc_graph, u_tab2)
                    # print(f"{tensor=}")
                    tc_elem = tensor.to_highest_weight()[0]
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
                    assert tc_elem.crystal_weight == tuple([a + b for a,b in zip_longest(dom.rc_graph.length_vector, u_tab_hw.length_vector, fillvalue=0)]), f"{tc_elem.crystal_weight=} vs {tuple([a + b for a,b in zip_longest(dom.rc_graph.length_vector, u_tab2.length_vector, fillvalue=0)])}"
                    high_weight_check = tuple([a for a, b in zip_longest(high_weight, tc_elem.crystal_weight, fillvalue=0)])
                    low_weight_check = tuple([a for a, b in zip_longest(rc_w.to_lowest_weight()[0].length_vector, tc_elem.crystal_weight, fillvalue=0)])
                    if tc_elem.crystal_weight == high_weight_check and tc_elem.to_lowest_weight()[0].crystal_weight == low_weight_check:
                        crystals[(rc_w, tc_elem)] = crystals.get(rc_w, 0) + 1
                        # print(f"{u=} {dom.perm=} {w=} {crystals=}")
                        highest_weights.add(tc_elem)
    return crystals

if __name__ == "__main__":

    ASx = FreeAlgebra(SchubertBasis)
    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)

    hw_tabs = set()
    for perm in perms:

        hw_tabs.update([RootTableau.from_rc_graph(rc.to_highest_weight()[0]) for rc in RCGraph.all_rc_graphs(perm, n - 1)])




    tot_suc = 0
    for w in perms:
        rc_w_coprods = {}
        good = False
        for hw_tab in hw_tabs:


            # if u.inv == 0:
            #     continue
            # needed?
            # if not u.bruhat_leq(dom.perm):
            #     continue

            # test compose tensor product, and homomorphism of RC graph ring


            coprod = ASx(w, n-1).coproduct()
            for ((d, _), (u, _)) in coprod:
                if d != hw_tab.perm:
                    continue
                crystals = decompose_tensor_product(hw_tab, u, n)

                rc_ring = RCGraphRing()

                tring = rc_ring @ rc_ring



                for (rc_w, tc_elem), coeff in crystals.items():
                    # if not rc_w.is:
                    #     continue
                    max_len = max(len(rc_w), len(tc_elem.factors[0]), len(tc_elem.factors[1]))
                    t_elem1, t_elem2 = tc_elem.factors
                    w_rc = rc_w.resize(max_len)
                    t_elem1 = t_elem1.resize(max_len)
                    t_elem2 = t_elem2.resize(max_len)

                    # mul_up = tring.one

                    # for index in range(max_len):
                    #     mul_up *= tring((RCGraph.one_row(len(t_elem1[index])), RCGraph.one_row(len(t_elem2[index]))))
                    #     # subtract off not there
                    #     for (rc1, rc2) in mul_up:
                    #         if rc1.rowrange(0, index + 1).perm != t_elem1.rowrange(0, index + 1).perm or rc2.rowrange(0, index + 1) != t_elem2.rowrange(0, index + 1).perm:
                    #             # rc_w_coprods[rc_w] = rc_w_coprods.get(rc_w, tring.zero) + coeff * tring((rc1, rc2))
                    #             mul_up -= tring((rc1, rc2))

                    # RAISE TO HIGHEST WEIGHT TO MAKE COASS?
                    # ???
                    # DON'T USE THE COEFF
                    # rc_w_coprods[w_rc] = rc_w_coprods.get(w_rc, tring.zero) + coeff * tring((t_elem1, t_elem2))
                    #if (t_elem1, t_elem2) not in rc_w_coprods.get(w_rc, tring.zero):
                    rc_w_coprods[w_rc] = rc_w_coprods.get(w_rc, tring.zero) + coeff * tring((t_elem1, t_elem2))
                    # if t_elem1.perm != t_elem2.perm and (t_elem2, t_elem1) not in rc_w_coprods.get(w_rc, tring.zero):
                    #     rc_w_coprods[w_rc] += tring((t_elem2, t_elem1))
                    # if not t_elem2.is_highest_weight:
                    #     rc_w_coprods[w_rc] += tring((t_elem2, t_elem1))
                #make_cocom = tring.zero
        # rc_w_coprods_old = {**rc_w_coprods}
        # rc_w_coprods = {}
        # for w_rc, val in rc_w_coprods_old.items():
        #     rc_w_coprods[w_rc] = rc_w_coprods.get(w_rc, tring.zero) + val
        #     for (rc1, rc2), coeff in val.items():
        #         assert coeff == 1, f"{coeff=} {rc1,rc2} {val=}"
        #         crystals = decompose_tensor_product(RootTableau.from_rc_graph(rc2.to_highest_weight()[0]), rc1.perm)
        #         for (rc_w, tc_elem), coeff in crystals.items():
        #         # if not rc_w.is:
        #         #     continue
        #             max_len = n - 1
        #             t_elem1, t_elem2 = tc_elem.factors
        #             w_rc2 = rc_w.resize(max_len)
        #             if w_rc2 != w_rc:
        #                 continue
        #             t_elem1 = t_elem1.resize(max_len)
        #             t_elem2 = t_elem2.resize(max_len)
        #             if (t_elem1, t_elem2) not in rc_w_coprods[w_rc]:
        #                 rc_w_coprods[w_rc] += tring((t_elem1, t_elem2))


        for rc, val in rc_w_coprods.items():
            if rc.perm != w:
                continue
            tens_ring = ASx@ASx
            check_elem = tens_ring.zero
            for (rc1, rc2), coeff in val.items():
                assert coeff == 1
                check_elem += tens_ring(((rc1.perm,len(rc1)), (rc2.perm,len(rc2))))
            diff = check_elem - coprod
            try:
                assert all(v == 0 for v in diff.values())
            except AssertionError:
                # print("A fail")
                # print(f"{diff=}")
                continue
            print(f"Coprod {rc.perm.trimcode}")
            pretty_print(rc)
            pretty_print(val)
            print("At least one success")
            good = True
            break
        assert good, f"COMPLETE FAIL {w=}"
        tot_suc += 1
                # prod = rc_ring.element_from_rc_graph(hw_tab.rc_graph) * rc_ring.element_from_rc_graph(tc_elem.factors[1])
                # for rc_g in prod:
                #     new_crystals[(rc_w, rc_g)] = new_crystals.get((rc_w, rc_g), 0) + coeff * prod[rc_g]

            # OK now we have the HW
            # product of rows equal to the hw of tensor product?

            # try:
            #     test_prod = Sx(u)*Sx(hw_tab.perm)
            #     for w in test_prod:
            #         any_cry = None
            #         for rcc in crystals:
            #             if rcc.perm == w:
            #                 any_cry = rcc
            #                 break
            #         coeff = crystals[any_cry]
            #         assert coeff == test_prod[w], f"Fail at coeff check, {u=} {w=} {test_prod[w]=} {coeff=}, {crystals=}"
            # except AssertionError as e:
            #     print(f"Failed coeff check at {u=} {w=} {dom.perm=}, expected {(Sx(dom.perm)*Sx(u)).get(w, 0)}, got {coeff}")
            #     print("Crystals:")
            #     for c in crystals:
            #         pretty_print(c)
                # print("Highest weights:")
                # for c in highest_weights:
                #     pretty_print(c)
            #     raise
            # print(f"Successful tensor decompose {dom.perm,u}")
    assert tot_suc == len(perms)
