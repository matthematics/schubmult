# # LR rule verification script

# import os
# import shutil
# import sys
# import time
# from functools import cache
# from itertools import zip_longest
# from json import dump, load
# from multiprocessing import Lock, Manager, Pool, Process, cpu_count

# from schubmult.rings.rc_graph_module import FA, RCGraph, try_lr_module

# all_rc_graphs = {}

# def generate_rc_graphs(perms, n):
#     global all_rc_graphs
#     for perm in perms:
#         all_rc_graphs[perm] = []
#         for length in range(n):
#             all_rc_graphs[perm].append({g: i for i, g in enumerate(tuple(sorted(RCGraph.all_rc_graphs(perm, length), key=lambda a: a.length_vector())))})

# def co_principal_perms(perm, length):
#     from schubmult import ASx, uncode

#     if len(perm.trimcode) == 0:
#         return []
#     lower_elem = ASx(uncode(perm.trimcode[1:]), length - 1)
#     upper_elem = ASx(uncode([perm.trimcode[0]]), 1) * lower_elem
#     return [key[0] for key in upper_elem.keys() if key[0] != perm]

# # @cache
# # def is_good(perm, rc1, rc2):
# #     from schubmult import FA, uncode
# #     if len(rc1) != len(rc2):
# #         return False
# #     if perm.inv == 0 and len(rc1) == 0:
# #         return True
# #     if uncode([a+b for a,b in zip(rc1.length_vector(), rc2.length_vector())]) != perm:
# #         return False
# #     if len(rc1) == 1:
# #         return True
# #     lower_perm = uncode(perm.trimcode[1:])
    
# #     if is_good(lower_perm, rc1_lower, rc2_lower) and (rc1, rc2) in (FA(perm.trimcode[0]).coproduct() * (rc1_lower @ rc2_lower)).value_dict:
# #         good = True
# #         break
# #     for perm2 in co_principal_perms(perm, len(rc1)):
# #         for rc01 in RCGraph.all_rc_graphs(rc1.perm, len(rc1)):
# #             for rc02 in RCGraph.all_rc_graphs(rc2.perm, len(rc2)):
# #                 if is_good(perm2, rc01, rc02) and (rc1 >= rc01 or rc2 >= rc02):
# #                     return False
# #     return True

# def main():
#     from schubmult import Permutation, uncode

#     # try:
    

#     n = int(sys.argv[1])
#     perms = Permutation.all_permutations(n)
#     #extra_perms = Permutation.all_permutations(n+n)
#     #generate_rc_graphs(extra_perms, n)
#     #perms.sort(key=lambda p: (len(p.trimcode), p.trimcode), reverse=True)
    
#     latest_seen_left = {perm: [{} for _ in range(n)] for perm in perms}
#     latest_seen_right = {perm: [{} for _ in range(n)] for perm in perms}
    
#     all_present_graphs = {}
#     for perm in perms:
#         # all_present_graphs[perm] = {}
#         modules = try_lr_module(perm)
#         # print(all_rc_graphs[perm])
#         # length = 0
#         graphs_to_check = FA(*perm.trimcode).coproduct() * (RCGraph()@RCGraph())
#         result = 0 * (RCGraph()@RCGraph())
#         for (rc1, rc2), coeff in graphs_to_check.items():
#             if rc1.perm.bruhat_leq(perm) and rc2.perm.bruhat_leq(perm):
#                 if is_good(perm, rc1, rc2):
#                     result += coeff * (rc1 @ rc2)

#         try:
#             assert all(v==0 for v in (modules - result).value_dict.values())
#         except AssertionError:
#             print(f"Failure for {perm.trimcode}")
#             print("Expected:")
#             print(modules)
#             print("Got:")
#             print(result)
#             continue
#         print(f"Success for {perm.trimcode}")
#         #for (rc1, rc2), coeff in modules.items():
#             # length = len(rc1)

#             # print(f"{perm.trimcode}: {rc1.perm.trimcode} {len(rc2.perm.trimcode)}")
#             # print(f"  {latest_seen_left[perm][len(rc1)].get(rc1.perm, -1)} * {latest_seen_right[perm][len(rc2)].get(rc2.perm, -1)} latest seen")

#             # index1 = all_rc_graphs[rc1.perm][len(rc1)][rc1]
#             # index2 = all_rc_graphs[rc2.perm][len(rc2)][rc2]
#             # #all_present_graphs[(rc1.perm, rc2.perm)].update([((rc1, rc2) for rc1 in all_rc_graphs(rc1.perm)[len(rc1)] if rc1 != rc1] ))
#             # print(f"{index1} * {index2} = {coeff}")
#             # #all_pairs = {(rc01, rc02) for rc1, rc2 in product(all_rc_graphs(rc1.perm)[len(rc1)], all_rc_graphs(rc2.perm)[len(rc2)]) if vector_sum(rc1.length_vector(), rc2.length_vector()) == vector_sum(rc01.length_vector(), rc02.length_vector())}

#             # #all_present_graphs[(rc1.perm, rc2.perm)].update()
#             # try:
#             #     assert not ((latest_seen_left[perm][len(rc1)].get(rc1.perm, 1000000), latest_seen_right[perm][len(rc2)].get(rc2.perm, 1000000000000000)) < (index1, index2))
#             # except AssertionError:
#             #     print(f"Order violation: {perm.trimcode}: {rc1.perm.trimcode} {len(rc2.perm.trimcode)}")
#             #     print("********VIOLATION*******")
#             #     print(f"  {latest_seen_left[perm][len(rc1)].get(rc1.perm)} * {latest_seen_right[perm][len(rc2)].get(rc2.perm)} latest seen")
#             #     print(f"  {index1} * {index2}")
#             #     print("********END VIOLATION*******")
#             # all_present_graphs[perm][(rc1.perm, rc2.perm)] = all_present_graphs[perm].get((rc1.perm, rc2.perm), set())
#             # all_present_graphs[perm][(rc1.perm, rc2.perm)].add((rc1, rc2))
#             # #print(f"  {index1} * {index2} = {coeff}")
#             # # latest_seen_left[perm][len(rc1)][rc1.perm] = index1
#             # # latest_seen_right[perm][len(rc2)][rc2.perm] = index2
#             # for co_perm in co_principal_perms(perm, length):
#             #     #print(f"    Co-principal: {co_perm.trimcode} {all_rc_graphs[co_perm][len(rc1)].get(rc1)}")
#             # #print(f"  {index1} * {index2} = {coeff}")
#             #     if len(co_perm) <= n:
#             #         assert co_perm.trimcode < perm.trimcode
#             #         latest_seen_left[co_perm][len(rc1)][rc1.perm] = index1# if (not latest_seen_left[co_perm][len(rc1)] or latest_seen_left[co_perm][len(rc1)][rc1.perm]) < index1 else latest_seen_left[co_perm][len(rc1)][rc1.perm]
#             #         latest_seen_right[co_perm][len(rc2)][rc2.perm] = index2# if (not latest_seen_right[co_perm][len(rc2)] or latest_seen_right[co_perm][len(rc2)][rc2.perm]) < index2 else latest_seen_right[co_perm][len(rc2)][rc2.perm]

#             # for (perm1, perm2), seen_set in all_present_graphs.items():
#             #     for rc1, index1 in all_rc_graphs[perm1][length].items():
#             #         for rc2, index2 in all_rc_graphs[perm2][length].items():
#             #             if uncode([a+b for a,b in zip_longest(rc1.length_vector(), rc2.length_vector(), fillvalue=0)]) == perm and (latest_seen_left[perm][len(rc1)].get(rc1.perm,-1),latest_seen_right[perm][len(rc2)].get(rc2.perm,-1)) < (index1, index2) and (rc1, rc2) not in seen_set:
#             #                 print(f"Missing graph for {perm.trimcode}: {rc1.perm.trimcode,rc1.length_vector()} {rc2.perm.trimcode,rc2.length_vector()} at length {length}")
#             # print("---------")
#     # for perm in perms:
#     #     length = len(perm.trimcode)
#     #     for co_perm in co_principal_perms(perm, length):
#     #         if len(co_perm) > n:
#     #             continue
#     #         bad = False
#     #         for (perm1, perm2), seen_set in all_present_graphs[co_perm].items():
#     #             for (rc01, rc02) in  seen_set:
#     #                 for rc1 in all_rc_graphs[rc01.perm][length]:
#     #                     for rc2 in all_rc_graphs[rc02.perm][length]:
#     #                         if uncode([a+b for a,b in zip(rc1.length_vector(), rc2.length_vector())]) == perm and not ((all_rc_graphs[rc1.perm][len(rc1)][rc01],all_rc_graphs[rc2.perm][len(rc2)][rc02]) > (all_rc_graphs[rc1.perm][len(rc1)][rc1], all_rc_graphs[rc2.perm][len(rc2)][rc2])) and (rc1, rc2) not in all_present_graphs[co_perm][(rc01.perm, rc02.perm)]:
#     #                             print(f"Missing graph for {co_perm.trimcode}: {rc1.perm.trimcode,rc1.length_vector()} {rc2.perm.trimcode,rc2.length_vector()} {(latest_seen_left[co_perm][len(rc1)].get(rc1.perm,-1),latest_seen_right[co_perm][len(rc2)].get(rc2.perm,-1))} at length {length}")
#     #                             bad = True
#     #         if bad:
#     #             print(f"Wah {co_perm.trimcode} is bad")


# if __name__ == "__main__":
#     main()
