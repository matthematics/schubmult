import copy
import itertools
import sys
from functools import cache
from itertools import zip_longest

import sympy
from sympy import pretty_print

from schubmult import CrystalGraphTensor, FreeAlgebra, NilPlactic, Permutation, RCGraph, RCGraphRing, RootTableau, SchubertBasis, Sx


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
    fyi = {}
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
                reduced_word = w_rc.reduced_word
                compatible_sequence = w_rc.compatible_sequence
                dom_tab = RootTableau.from_rc_graph(dom)
                w_tab = RootTableau.from_rc_graph(w_rc)
                u = u_rc.perm
                if tensor_hw.crystal_weight == high_weight and tensor_lw.crystal_weight == low_weight:# and (u_rc2.perm.minimal_dominant_above() == u_rc2.perm or w_rc.perm.minimal_dominant_above() != w_rc.perm):
                    found_one = False
                    for subword in all_reduced_subwords(w_rc.reduced_word, u):
                        print("w_tab")
                        pretty_print(w_tab)
                        print(f"{w_rc.reduced_word=}")
                        pretty_print(dom_tab)
                        roots = [w_tab.perm.right_root_at(index, word=w_rc.reduced_word) for index in subword]
                        grid = copy.deepcopy(w_tab._root_grid)
                        for box in w_tab.iter_boxes:
                            print(box)
                            if grid[box][0] in roots:
                                grid[box] = (grid[box][0], MarkedInteger(grid[box][1]))
                        u_tab = RootTableau(grid)
                        # compatible_seq = [MarkedInteger(a) if index in subword else a for index, a in enumerate(compatible_sequence)]
                        # u_tab = RootTableau.root_insert_rsk(reduced_word, compatible_seq)
                        # skip = False
                        # for box in dom_tab.iter_boxes:
                        #     if isinstance(u_tab[box][1], MarkedInteger):
                        #         skip = True
                        #         break
                        # if skip:
                        #     continue
                        last_inv = 1000

                        d_tab = RootTableau.from_rc_graph(w_rc)
                        froff = False
                        while u_tab.perm.inv < last_inv:
                            last_inv = u_tab.perm.inv
                            for box in u_tab.iter_boxes_row_word_order:
                                if not isinstance(u_tab[box][1], MarkedInteger):
                                    u_tab_test = u_tab.delete_box(box)
                                    if u_tab_test is not None:
                                        u_tab = u_tab_test
                                        break
                                # else:
                                #     try:

                                #         d_tab_test = d_tab.up_jdt_slide(*box, force=True)
                                #         if d_tab_test is not None:
                                #             d_tab = d_tab_test
                                #     except Exception:
                                #         froff = False
                                #         print("Couldn't up jdt")
                                #         pretty_print(d_tab)
                                #         print(f"{box=}")
                        if u_tab.perm.inv > u_rc.perm.inv:
                            # didn't make it
                            print("No make")
                            print(u_tab)
                            continue

                        # u_tab = u_tab.anti_rectify(w_tab.rows, w_tab.cols)
                        # pretty_print(u_tab)
                        # print("ANTIEAF")
                        # assert u_tab.rows >= w_tab.rows
                        # assert u_tab.cols >= w_tab.cols
                        # for box in dom_tab.iter_boxes:
                        #     if u_tab[box] is not None:
                        #         print("THESTINK")
                        #         d_tab = None
                        #         break
                        # if d_tab.rc_graph.resize(n - 1) == dom:
                        #     print("HOORAY")
                        #     found_one = True

                        #d_tab = d_tab.anti_rectify()
                        # d_tab = RootTableau.from_rc_graph(w_rc.to_highest_weight()[0])
                        # last_inv = 1000
                        # while d_tab.perm.inv < last_inv:
                        #     last_inv = d_tab.perm.inv
                        #     for box in d_tab.iter_boxes:
                        #         if dom_tab[box] is None:
                        #             try:
                        #                 d_tab_test = d_tab.delete_box(box)
                        #                 if d_tab_test is not None:
                        #                     d_tab = d_tab_test
                        #             except Exception:
                        #                 print("Bleh")
                        # print(root_grid)
                        # try:
                        #     d_tab = RootTableau(root_grid)
                        #     d_tab = d_tab.anti_rectify()
                        # except Exception:
                        #     print("oops")
                        #     d_tab = w_tab
                        # # unrectify
                        # d_tab2 = RootTableau.from_rc_graph(dom)
                        # while len(list(d_tab2.iter_outer_corners)) > 0:
                        #     for box in d_tab2.iter_outer_corners:
                        #         d_tab2 = d_tab2.up_jdt_slide(*box)
                        #         break
                        # new_d_tab_grid = copy.deepcopy(d_tab._root_grid)
                        # for box2 in d_tab2.iter_boxes:
                        #     if new_d_tab_grid[box2] is None:
                        #         print("Could not")
                        #         pretty_print(d_tab)
                        #         break
                        #     new_d_tab_grid[box2] = (d_tab2[box2][0],d_tab2[box2][1])
                        # try:
                        #     d_tab = RootTableau(new_d_tab_grid).rectify()
                        # except Exception:
                        #     print("Couldn't do it")

                        print("u_tab")
                        pretty_print(u_tab)
                        u_hw_rc = u_tab.rc_graph.resize(n-1)
                        if u_hw_rc.perm != u:
                            print("Got")
                            pretty_print(u_hw_rc)
                            print("Perm is not")
                            print(u)
                            print(f"{u.antiperm=}")
                            continue
                        if u_hw_rc == u_rc:
                            print("Got the identical crystal")
                            # if d_tab.perm.inv == dom.perm.inv:
                            #     print("THEPERMIVN")
                            #     pretty_print(d_tab)
                            #     assert d_tab.rc_graph.resize(n-1).length_vector == dom.length_vector
                            #     print("EHHOG")
                            crystals[w_rc] = crystals.get(w_rc, set())
                            crystals[w_rc].add(tensor)
                            fyi[w_rc] = fyi.get(w_rc, set())
                            fyi[w_rc].add((tensor, u_tab))
                            # d_tab_set[w_rc] = d_tab_set.get(w_rc, set())
                                # d_tab_set[w_rc].add(d_tab)
                    # if found_one:
                    #     print("FOOOOOFUND")
                    #     fyi[w_rc] = fyi.get(w_rc, set())
                    #     fyi[w_rc].add((tensor, found_one))
                # for subword in all_reduced_subwords(w_rc.reduced_word, dom.perm):

                #         w_rc2 = w_rc#.to_highest_weight()[0]
                #         print(f"{w_rc2.reduced_word=}")
                #         pretty_print(dom_tab)
                #         print("w_tab2 DOMDOM")
                #         w_tab2 = RootTableau.from_rc_graph(w_rc2)
                #         pretty_print(w_tab2)
                #         roots = [w_tab.perm.right_root_at(index, word=w_rc2.reduced_word) for index in subword]
                #         grid = copy.deepcopy(w_tab2._root_grid)
                #         for box in w_tab2.iter_boxes:
                #             print(box)
                #             if grid[box][0] in roots:
                #                 grid[box] = (grid[box][0], MarkedInteger(grid[box][1]))
                #         d_tab = RootTableau(grid)
                #         # compatible_seq = [MarkedInteger(a) if index in subword else a for index, a in enumerate(compatible_sequence)]
                #         # u_tab = RootTableau.root_insert_rsk(reduced_word, compatible_seq)
                #         # skip = False
                #         # for box in dom_tab.iter_boxes:
                #         #     if isinstance(u_tab[box][1], MarkedInteger):
                #         #         skip = True
                #         #         break
                #         # if skip:
                #         #     continue
                #         last_inv = 1000

                #         #d_tab = RootTableau.from_rc_graph(w_rc)
                #         while d_tab.perm.inv < last_inv:
                #             last_inv = d_tab.perm.inv
                #             for box in d_tab.iter_boxes:
                #                 if not isinstance(d_tab[box][1], MarkedInteger):
                #                     d_tab_test = d_tab.delete_box(box)
                #                     if d_tab_test is not None:
                #                         d_tab = d_tab_test
                #                         break

                #         if d_tab.perm.inv > dom.perm.inv:
                #             # didn't make it
                #             print("No make")
                #             print(d_tab)
                #             continue

                #         d_tab = d_tab.rectify()
                #         # if d_tab.rc_graph.resize(n - 1) == dom:
                #         #     print("HOORAY")
                #         #     found_one = True

                #         #d_tab = d_tab.anti_rectify()
                #         # d_tab = RootTableau.from_rc_graph(w_rc.to_highest_weight()[0])
                #         # last_inv = 1000
                #         # while d_tab.perm.inv < last_inv:
                #         #     last_inv = d_tab.perm.inv
                #         #     for box in d_tab.iter_boxes:
                #         #         if dom_tab[box] is None:
                #         #             try:
                #         #                 d_tab_test = d_tab.delete_box(box)
                #         #                 if d_tab_test is not None:
                #         #                     d_tab = d_tab_test
                #         #             except Exception:
                #         #                 print("Bleh")
                #         # print(root_grid)
                #         # try:
                #         #     d_tab = RootTableau(root_grid)
                #         #     d_tab = d_tab.anti_rectify()
                #         # except Exception:
                #         #     print("oops")
                #         #     d_tab = w_tab
                #         # # unrectify
                #         # d_tab2 = RootTableau.from_rc_graph(dom)
                #         # while len(list(d_tab2.iter_outer_corners)) > 0:
                #         #     for box in d_tab2.iter_outer_corners:
                #         #         d_tab2 = d_tab2.up_jdt_slide(*box)
                #         #         break
                #         # new_d_tab_grid = copy.deepcopy(d_tab._root_grid)
                #         # for box2 in d_tab2.iter_boxes:
                #         #     if new_d_tab_grid[box2] is None:
                #         #         print("Could not")
                #         #         pretty_print(d_tab)
                #         #         break
                #         #     new_d_tab_grid[box2] = (d_tab2[box2][0],d_tab2[box2][1])
                #         # try:
                #         #     d_tab = RootTableau(new_d_tab_grid).rectify()
                #         # except Exception:
                #         #     print("Couldn't do it")
                #         print("d_tab")
                #         pretty_print(d_tab)
                #         dom_hw_rc = d_tab.rc_graph.resize(n-1)
                #         if dom_hw_rc.perm != dom.perm:
                #             print("Got")
                #             pretty_print(dom_hw_rc)
                #             print("Perm is not")
                #             print(u)
                #             print(f"{dom.perm.antiperm=}")
                #             continue
                #         if dom_hw_rc == dom:
                #             print("Got the identical crystal DOMDOM")
                #             fyi[w_rc] = fyi.get(w_rc, set())
                #             fyi[w_rc].add((w_rc, d_tab))
                #             # d_tab_set[w_rc] = d_tab_set.get(w_rc, set())
                #                 # d_tab_set[w_rc].add(d_tab)

    try:

        assert len(crystals) == 1
    except AssertionError:
        print("Error: More than one crystal found FOR ")
        pretty_print(u_rc)
        for rc, tensor_set in fyi.items():
            pretty_print(rc)
            for tensor, d_tab in tensor_set:
                if isinstance(tensor, RCGraph):
                    print("THISISRC")
                pretty_print(tensor)
                print("d_tab - maybe can infuse?")
                pretty_print(d_tab)
                print("rc for d_tab")
                if isinstance(d_tab, RootTableau):
                    pretty_print(d_tab.rc_graph)
                print("dom")
                pretty_print(dom)
                # print("high weight")
                # pretty_print(rc.to_highest_weight()[0])
                # pretty_print(tensor.to_highest_weight()[0])
                # print("low weight")
                # pretty_print(rc.to_lowest_weight()[0])
                # pretty_print(tensor.to_lowest_weight()[0])

    return crystals

def is_decomposable(w):
    for i in range(1, len(w) - 1):
        coset, w_J = w.coset_decomp(*list(range(1, i + 1)),*list(range(i + 2, len(w))))
        if coset.inv == 0 and set(w_J.code[:i+1]) != {0} and set(w_J.code[i+2:]) != {0}:
            return True
    return False


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("n",action="store",type=int)
    parser.add_argument("--case",action="append",type=int,nargs="+")
    parser.add_argument("-",nargs="+",action="append",type=int,dest="case")
    args = parser.parse_args()
    ASx = FreeAlgebra(SchubertBasis)
    n = args.n

    case = args.case
    # print(case)
    # k = int(sys.argv[2])
    # assert n - 1 == k

    perms = Permutation.all_permutations(n)

    hw_tabs = set()
    for perm in perms:
        if perm.minimal_dominant_above() != perm:# or perm != Permutation.ref_product(3,2,1,2):
            continue
        if case is not None:
            if len(case) == 1 and perm != Permutation.ref_product(*case):
                continue
            elif len(case) == 2 and perm != Permutation.ref_product(*case[0]):
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
        if case is not None and len(case) > 1 and Permutation.ref_product(*case[1]) != v:
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
                if the_rc.is_principal:
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
            #the_sum += coeff * rc.polyvalue(Sx.genset)
            if rc.is_principal:
                the_sum2 += coeff * Sx(rc.perm)
            #the_sum += coeff * Sx(rc.perm)
        prod2 = Sx(the_sum)
        print(f"{prod=}")
        print(f"{prod2=}")
        print(f"{the_sum2=}")
        #assert prod2 - prod == Sx.zero
        assert prod == the_sum2
        print(prod)


