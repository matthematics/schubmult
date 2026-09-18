"""RC graph ring comodule"""

from schubmult import *

if __name__ == "__main__":
    r = RCGraphRing()

    

    def cope(rc_graph):
        if len(rc_graph) == 0:
            return r.one @ FA()

        if len(rc_graph) == 1:
            result = 0
            for p in range(rc_graph.perm.inv + 1):
                result += r.monomial(p) @ FA(rc_graph.perm.inv - p)
            return result
        lower_result = cope(rc_graph.rowrange(1))
        result = 0
        for (rc, word), coeff in lower_result.items():
            for p in range(len(rc_graph[0]) + 1):
                # result += r.monomial(p) @ FA(word + (rc_graph.perm.inv - p,))
                up_rcs = r.monomial(p)* r(rc)
                for up_rc, coeff2 in up_rcs.items():
                    pants_fingers = FA(len(rc_graph[0]) - p, *word).change_basis(SchubertBasis)
                    for (perm0, length), coeff3 in pants_fingers.items():
                        if (Sx(up_rc.perm) * Sx(perm0)).get(rc_graph.perm, 0) != 0 and len(RCGraph.all_rc_graphs(perm0, length, weight=tuple([len(rc_graph[0]) - p, *tuple([b - a for b, a in zip(rc_graph.length_vector[1:],word)])]))) != 0:
                            result += r(up_rc) @ FA(len(rc_graph[0]) - p, *word)
                            
        return result

    def cope_cheat(rc_graph):
        cheat_elem = FA(*rc_graph.length_vector)
        the_prodprod = cheat_elem.coproduct()
        ret = r.zero @ FA()
        for (word1, word2), coeff in the_prodprod.items():
            ret += coeff * r.monomial(*word1) @ FA(*word2)
        return ret

    def cope_asx(rc_graph):
        if len(rc_graph) == 0:
            return r.one @ ASx([])

        if len(rc_graph) == 1:
            result = 0
            for p in range(rc_graph.perm.inv + 1):
                result += r.monomial(p) @ ASx(uncode([rc_graph.perm.inv - p]), 1)
            return result
        result = 0
        banging_coprod = ASx(rc_graph.perm, len(rc_graph)).coproduct()
        for ((perm1, length1), (perm2, length2)), coeff in banging_coprod.items():
            potato = False
            for rc in RCGraph.all_rc_graphs(perm1, length1):
                
                if all([rc_graph.length_vector[i] >= rc.length_vector[i] for i in range(len(rc_graph))]):
                    diff_vector = tuple([rc_graph.length_vector[i] - rc.length_vector[i] for i in range(len(rc_graph))])
                    for rc2 in RCGraph.all_rc_graphs(perm2, length2, weight=diff_vector):
                        if all([rc_graph.length_vector[i] >= rc2.length_vector[i] for i in range(len(rc_graph))]):
                            result += r(rc) @ ASx(perm2, length2)
                            potato = True
                            
                if potato:
                    break
                
            
                
        # lower_result = cope_asx(rc_graph.rowrange(1))
        # result = 0
        # for (rc, (perm, length)), coeff in lower_result.items():
        #     for p in range(len(rc_graph[0]) + 1):
        #         # result += r.monomial(p) @ FA(word + (rc_graph.perm.inv - p,))
        #         up_rcs = r.monomial(p)* r(rc)
        #         for up_rc, coeff2 in up_rcs.items():
        #             pants_fingers = FA(*word).change_basis(SchubertBasis)
        #             for (perm0, length), coeff3 in pants_fingers.items():
        #                 if (Sx(up_rc.perm) * Sx(perm0)).get(rc_graph.perm, 0) != 0 and len(RCGraph.all_rc_graphs(perm0, len(word), weight=word)) != 0:
        #                     result += coeff * coeff2 * r(up_rc) @ FA(len(rc_graph[0]) - p, *word)
                            
        return result

    def full_cope(rc_elem):
        result = 0
        for rc, coeff in rc_elem.items():
            result += coeff * cope_cheat(rc)
        return result

    import sys
    from sympy import pretty_print

    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)
    for perm in perms:
        pants = 0
        for rc in RCGraph.all_rc_graphs(perm):
            #pretty_print(rc)
            cop = full_cope(r(rc))
            #pretty_print(cop)

            # test comodule condition
            comod_elem1 = 0
            for (rc1, word1), coeff1 in cop.items():
                #comod_elem1 += coeff1 * r(rc1)@ (FA(*word1).coproduct())
                comod_elem1 += coeff1 * r.monomial(*rc1.length_vector)@ (FA(*word1).coproduct())

            comod_elem2 = 0
            for (rc1, word1), coeff1 in cop.items():
                comod_elem2 += coeff1 * full_cope(r(rc1)) @ FA(*word1)
            try:
                assert comod_elem1.almosteq(comod_elem2), f"Comodule condition failed for {rc} with coproduct {cop}\n{comod_elem1-comod_elem2}"
            except AssertionError as e:
                # print(e)
                # pretty_print(comod_elem1)
                # pretty_print(comod_elem2)
                pants += 1
        print(f"Checked comodule condition for all RC graphs of {perm}, failed {pants} times out of {len(RCGraph.all_rc_graphs(perm))}")