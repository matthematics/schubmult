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
@cache
def decompose_tensor_product(dom, u_rc, length, n):
    # global hw_rc_sets
    u = u_rc.perm
    u_rc = u_rc.resize(n-1)
    crystals = {}
    highest_weights = set()
    perm_set = set((Sx(u)*Sx(dom.perm)).keys())
    for w in perm_set:
        # if not u.bruhat_leq(w):
        #     continue
        # if not dom.perm.bruhat_leq(w):
        #     continue

        # print(f"Moving on to {u=} {w=} {dom.perm=}")
        if w not in hw_rc_sets:
            hw_rc_sets[w] = set()
            for rc_w in RCGraph.all_rc_graphs(w, n - 1):
                # pretty_print(rc_w)
                hw_rc_sets[w].add(rc_w.to_highest_weight(length=length)[0])
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
                if u_hw_rc != u_rc:
                    continue
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
                        crystals[(rc_w, tc_elem)] = crystals.get(rc_w, 0) + 1
                        # print(f"{u=} {dom.perm=} {w=} {crystals=}")
                        highest_weights.add(tc_elem)
    return crystals

if __name__ == "__main__":

    ASx = FreeAlgebra(SchubertBasis)
    n = int(sys.argv[1])
    k = int(sys.argv[2])

    perms = Permutation.all_permutations(n)

    hw_tabs = set()
    for perm in perms:

        hw_tabs.update([RootTableau.from_rc_graph(rc.to_highest_weight(length=k)[0]) for rc in RCGraph.all_rc_graphs(perm, n - 1)])



    rc_ring = RCGraphRing()
    tot_suc = 0
    # u times v
    print("NOTE THIS IS CORRECT AND AN ASSOCIATIVE ACTION FOR DOMINANT PERMS")
    the_schubs = {}
    # top_k_dom = {}
    used = {}
    for hw_tab0, v in itertools.product(hw_tabs, perms):
        if hw_tab0.perm.inv == 0 or len(v.descents()) > 1 or k - 1 not in v.descents() or not set(v.trimcode).issubset({0,1}):
            continue
    #     # rc_w_coprods = {}
    #     # good = False
    #     rc = hw_tab0.rc_graph.resize(n-1)
    #     div_perm = Permutation([])
    #     if rc.vertical_cut(k)[0].perm.inv < Permutation.w0(k + 1).inv or rc.vertical_cut(k)[0].perm != rc.vertical_cut(k)[0].perm.minimal_dominant_above():
    #         the_perm = rc.perm
    #         found_any = True
    #         while found_any:
    #             found_any = False
    #             for i in range(k - 1):
    #                 if the_perm[i] < the_perm[i + 1]:
    #                     the_perm = the_perm.swap(i, i + 1)
    #                     div_perm = div_perm.swap(i, i + 1)
    #                     found_any = True
    #                     break
    #         new_rc = None
    #         for rc2 in RCGraph.all_rc_graphs(the_perm, n-1):
    #             if rc2.rowrange(k) == rc.rowrange(k):
    #                 new_rc = rc2
    #                 break
    #         assert new_rc is not None
    #         top_k_dom[new_rc] = top_k_dom.get(new_rc, set())
    #         if rc.perm in top_k_dom[new_rc]:
    #             continue
    #         top_k_dom[new_rc].add(rc.perm)
    #         hw_tab = RootTableau.from_rc_graph(new_rc)
    #     else:
    #         hw_tab = hw_tab0
        hw_tab = hw_tab0
        print("hw_tab")

        pretty_print(hw_tab.rc_graph)
        print(f"Da cut at {k}")
        the_cut0, the_cut1 = hw_tab.rc_graph.resize(n-1).vertical_cut(k)

        div_perm = (~the_cut0.perm)*the_cut0.perm.minimal_dominant_above()
        min_dom_graph = RCGraph.principal_rc(the_cut0.perm.minimal_dominant_above(), n-1)

        # find an exchange property path from min_dom_graph to the crystal of the_cut0
        def is_subgraph(rc1, rc2):
            for i in range(len(rc1)):
                if len(rc2[i]) < len(rc1[i]):
                    return False
                for j in range(len(rc1[i])):
                    if rc1[i][j] not in rc2[i]:
                        return False
            return True
        exchg_seq = []
        g = min_dom_graph.resize(len(the_cut0))
        while g != the_cut0:
            for d in sorted(g.perm.descents(), reverse=True):
                g0 = g.exchange_property(d + 1)
                if is_subgraph(g0, the_cut0):
                    exchg_seq.append(d + 1)
                    g = g0.to_lowest_weight(length = k)[0]
                    break

        # used[(min_dom_graph, the_cut1)] = used.get((min_dom_graph,the_cut1), set())
        # if hw_tab.perm in used[(min_dom_graph, the_cut1)]:
        #     continue
        # used[(min_dom_graph, the_cut1)].add(hw_tab.perm)
        # parallel exchange property
        for rc_v in RCGraph.all_rc_graphs(v, n - 1):
            crystals = decompose_tensor_product(RootTableau.from_rc_graph(min_dom_graph), rc_v, length=k, n=n)

            print("Product:")
            pretty_print(hw_tab0.rc_graph)
            print("and")
            pretty_print(rc_v)
            # MUST MODIFY SM. RULE: WEIGHT PRESERVING, DIVDIFF from div_perm
            # THIS IS CRYSTAL LEVEL

            # sm = rc_ring.from_dict({(k[0]: v for k, v in crystals.items()})
            sm = rc_ring.zero
            for (the_rc, tc_elem), coeff in crystals.items():
                # ALMOST CORRECT BUT WE HAVE SOME TWOS
                assert coeff == 1
                trim_down = the_rc
                bad = False
                for d in exchg_seq:
                    if d - 1 in trim_down.perm.descents():
                        trim_down = trim_down.exchange_property(d)
                    else:
                        bad = True
                        break
                assert not bad
                sm += rc_ring(trim_down)
                # permo = the_rc.perm
                # assert len(permo.trimcode) <= k
                # if (permo * (~div_perm)).inv != permo.inv - div_perm.inv:
                #     continue
                # new_perm = permo * (~div_perm)
                # tried = set()
                # if new_perm not in hw_rc_sets:
                #     hw_rc_sets[new_perm] = set()
                #     for rc_w in RCGraph.all_rc_graphs(new_perm, n - 1):
                #         # pretty_print(rc_w)
                #         hw_rc_sets[new_perm].add(rc_w.to_highest_weight(length=k)[0])
                # for rc_new in hw_rc_sets[new_perm]:
                #     # if rc_new.rowrange(k) == the_rc.rowrange(k):
                #     #     sm += rc_ring.from_dict({rc_new: 1})
                #     actual_rc_new_set = rc_new.vertical_cut(k)[0].prod_with_rc(the_cut1)
                #     true_set = {rc0.to_highest_weight(length=k)[0] for rc0 in actual_rc_new_set}
                #     for rc_add in true_set:

                #         if rc_add.perm in (Sx(hw_tab0.perm) * Sx(v)) and rc_add not in sm:
                #             sm += rc_ring.from_dict({rc_add: 1})
            pretty_print(sm)
            the_schubs[(hw_tab0.perm, v)] = the_schubs.get((hw_tab0.perm, v), rc_ring.zero) + sm
            break
    Permutation.print_as_code = True
    for (u, v), val in the_schubs.items():
        print(f"{u} * {v}=")
        pretty_print(val)
        prod = Sx(u) * Sx(v)
        for rc, coeff in val.items():
            assert prod[rc.perm] == coeff
        print(prod)


