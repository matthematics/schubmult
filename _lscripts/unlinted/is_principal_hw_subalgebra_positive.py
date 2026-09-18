from schubmult import *
from schubmult.rings.combinatorial import *
from schubmult.rings.free_algebra import *
from sympy import pretty_print

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    r = RCGraphRing()
    for perm in Permutation.all_permutations(n):
        if len(perm.trimcode) <= 1:
            continue
        #for length in range(max(2,len(perm.trimcode)), n):
        length = len(perm.trimcode)
        # if length <= 1:
        #     continue
        #for rc in RCGraph.all_hw_rcs(perm, length):
        for rc in {RCGraph.principal_rc(perm, length)}:
            #rc = rc.to_highest_weight()[0]
            #elem = r.monomial(*rc.length_vector)  # ASx(perm, length).change_basis(WordBasis)
            elem = r.from_free_algebra_element(ASx(rc.perm, len(rc)))

            # for i in range(1, length):
            #     #elem -= elem.apply_to_keys(lambda k: k.chute_lower(i)).apply_to_keys(lambda k: k.chute_raise(i))
            #     elem -= elem.lowering_operator(i).raising_operator(i)
            #     if not rc.is_highest_weight:
            #         rc0 = rc.raising_operator(i)
            #         if rc0 is not None:
            #             elem = elem.raising_operator(i).lowering_operator(i)
            # elem_last = r.zero
            # rc1 = rc
            # found = True
            # while found :
            #     #elem_last = elem
            #     found = False
            #     for i in range(1, length):
            #         #elem -= elem.apply_to_keys(lambda k: k.chute_lower(i)).apply_to_keys(lambda k: k.chute_raise(i))
            #         rc0 = rc1.chute_raise(i)
            #         if rc0 is not None:
            #             rc1 = rc0
            #             found = True
            #             elem = elem.apply_to_keys(lambda k: k.chute_raise(i))#.apply_to_keys(lambda k: k.chute_lower(i))
            #             break
                    #elem = elem.lowering_operator(i).raising_operator(i)
            #elem_last = r.zero
            
            # found = True
            # rc_move = rc
            # while found:
            #     found = False
            #     for i in range(1, length):
            #         rc0 = rc_move.chute_raise(i)
            #         if rc0 is not None:
            #             rc_move = rc0
            #             found = True
            #             elem = elem.apply_to_keys(lambda k: k.chute_raise(i))#.apply_to_keys(lambda k: k.chute_lower(i))
            #             break

            # assert rc_move.is_highest_weight, f"Error: raising operators applied to {perm.trimcode} did not yield a highest weight RC graph, got {rc_move}"
            # assert rc_move in elem.keys(), f"Error: raising operators applied to {perm.trimcode} did not yield an RC graph in the resulting element, got {rc_move} with element {elem}"

            found = True
            while found :
                #elem_last = elem
                found = False
                for rc2 in list(elem.keys()):
                    for i in range(1, length):
                        #the_next = rc2.chute_lower(i)
                        the_next = rc2.lowering_operator(i)
                        if the_next is not None:
                            spitnb = elem.apply_to_keys(lambda k: k.lowering_operator(i)).apply_to_keys(lambda k: k.raising_operator(i))
                            assert spitnb != r.zero, f"Error: lowering and raising operators applied to {perm.trimcode} yielded zero, expected a nonzero element corresponding to {the_next}"
                            elem -= spitnb
                            found = True
                            assert rc2 not in elem.keys() or elem[rc2] == 0, f"Error: raising operators applied to {perm.trimcode} did not eliminate the non-highest weight RC graph {rc2}, got element {elem}"
                            break
                    if found:
                        break
                        # rc0 = rc1.chute_lower(i)
                        # if rc0 is not None:
                        #     rc1 = rc0
                        #     found = True
                        #     elem = elem.apply_to_keys(lambda k: k.chute_lower(i))#.apply_to_keys(lambda k: k.chute_raise(i))
                        #     break
            # while not elem_last.almosteq(elem):
            #     elem_last = elem
            #     for i in range(1, length):
            #         #elem -= elem.apply_to_keys(lambda k: k.chute_lower(i)).apply_to_keys(lambda k: k.chute_raise(i))
            #         rc0 = rc1.chute_lower(i)
            #         if rc0 is not None:
            #             rc1 = rc0
            #             elem = elem.apply_to_keys(lambda k: k.chute_lower(i))#.apply_to_keys(lambda k: k.chute_raise(i))
            #         #elem = elem.lowering_operator(i).raising_operator(i)
            assert len(elem) == 1, f"Error: raising operators applied to {perm.trimcode} did not yield a single RC graph, got {elem}"
            
        # rc = RCGraph.principal_rc(perm)#.to_highest_weight()[0]
        # length = len(rc)
        
        # found = True
        # while found:
        #     found = False
        #     for i in range(1, len(rc)):
        #         t_elem = elem.lowering_operator(i).raising_operator(i)
        #         if not t_elem.almosteq(elem) and not t_elem.almosteq(r.zero):
        #             elem -= t_elem
        #             found = True
        #     pretty_print(elem)
                    
        #     # rc_start = rc
        #     # stack = [rc_start]
        #     # while stack:
        #     #     the_rc = stack.pop()
        #     #     for i in range(1, len(rc)):
        #     #         test_rc = the_rc.raising_operator(i)
        #     #         if test_rc is None:                            
        #     #             elem -= elem.raising_operator(i)
        #     #         else:
        #     #             stack.append(test_rc)
        # #assert next(iter(elem.keys())) == rc_hw, f"Error: raising operators applied to {perm.trimcode} did not yield the expected RC graph {rc}, got {elem}"
        #assert len(elem.keys()) == 1, f"Error: raising operators applied to {perm.trimcode} did not yield a single RC graph, got {elem}"
        #print("BAinging")
            #elem = ASx(perm, length).change_basis(WordBasis)
        #length = len(rc)
            # elem = r.from_free_algebra_element(ASx(perm, length))
            # rc_hw, raise_seq = rc.to_highest_weight()
            # for i in raise_seq:
            #     elem = elem.raising_operator(i)
            # # for i in range(1, length):
            # #     sub_elem = elem
            # #     n = 0
            # #     while True:
            # #         elem0 = sub_elem.raising_operator(i)
            # #         if elem0 == r.zero:
            # #             break
            # #         n += 1
            # #         sub_elem = elem0
            # #         for i in range(n):
            # #             elem0 = elem0.lowering_operator(i)                
            # #         elem -= elem0
            # elem = elem.reverse_raise_seq(raise_seq)
            # # while not all(rc0.is_highest_weight for rc0,v in elem.items() if v != 0):
            # #     for i in range(1, length):
            # #         elem -= elem.raising_operator(i)
            #         # if elem != r.zero:
            #         #     elem0 = elem.lowering_operator(i)
            #         #     if elem0 == r.zero:
            #         #         break
            #         #     elem = elem0
            # assert elem != r.zero, f"Error: raising operators applied to {perm.trimcode} yielded zero, expected a highest weight element corresponding to {rc}"
            # assert all(rc0.is_principal for rc0,v in elem.items() if v != 0), f"Error: raising operators applied to {perm.trimcode} did not yield only highest weight RC graphs, got {elem}"
            # pretty_print(elem)
            # # lw_
            # # rc_sum = 0
            # # for word, coeff in elem.items():
            # #     rc_sum += coeff * r.monomial(*word)
            # # #rc_sum_hw = sum(coeff * r(rc.to_highest_weight()[0]) for rc, coeff in rc_sum.items())
            # # found = True
            # # while found:
            # #     found = False
            # #     for i in range(1, length):
            # #         t_elem = rc_sum.lowering_operator(i).raising_operator(i)
            # #         if not t_elem.almosteq(rc_sum) and not t_elem.almosteq(r.zero):
            # #             rc_sum -= t_elem
            # #             found = True
            # #     #pretty_print(rc_sum)
            # # assert len(rc_sum) == 1, f"Expected exactly one term in the coproduct of {perm.trimcode}, got {rc_sum}"
            # # assert all(v >= 0 for _, v in rc_sum.items()), f"Negative coefficient in coproduct of {perm.trimcode}, got {rc_sum}"
        print("bingo")