from schubmult import Permutation, elem_sym_perms, elem_sym_perms_q, q_vector, phi_d, tau_d, RCGraph, SchubertBasis, ASx, WordBasis, FA, Sx, RCGraphRing, Plactic, RootTableau
from schubmult.abc import x
from symengine import S

def antipode(fa_elem):
    sc_elem = fa_elem.change_basis(WordBasis)
    anti_elem = {}
    for sq, coeff in sc_elem.items():
        anti_elem[tuple(reversed(sq))] = (S.NegativeOne ** (len(sq))) * coeff
    return FA.from_dict(anti_elem).change_basis(fa_elem.ring._basis)

# def reverse_kk_with_rows_and_reflections(rc: RCGraph, d, row_reflection_pairs):
#     """
#     Reverse KK insertion using pairs of (row, reflection).
    
#     Key insight: Each (target_row, (a,b)) means:
#     - There's an element with root (a,b) at or below target_row
#     - It was pushed down during insertion and should be moved back UP to rows < target_row
    
#     row_reflection_pairs: list of (row, (a, b)) tuples from kogan_kumar_insert
#     """
#     if len(row_reflection_pairs) == 0:
#         return rc
    
#     working_rc = rc
    
#     # Process from largest target_row to smallest (bottom to top)
#     sorted_pairs = sorted(row_reflection_pairs, key=lambda x: x[0] if x[0] is not None else 0, reverse=True)
    
#     for target_row, (a, b) in sorted_pairs:
#         if target_row is None:
#             continue
            
#         # Find the element with root (a,b) at or below target_row
#         found_row, found_col = None, None
#         for row in range(target_row, len(working_rc) + 1):
#             for col in range(1, working_rc.cols + 20):
#                 if working_rc.has_element(row, col):
#                     q, s = working_rc.right_root_at(row, col)
#                     if (q, s) == (a, b):
#                         found_row, found_col = row, col
#                         break
#             if found_row is not None:
#                 break
        
#         if found_row is None:
#             continue  # Element not found, skip
            
#         # Remove element from current position
#         working_rc = working_rc.toggle_ref_at(found_row, found_col)
        
#         # Find where it should go in rows < target_row
#         for dest_row in range(1, target_row):
#             for dest_col in range(1, working_rc.cols + 20):
#                 if not working_rc.has_element(dest_row, dest_col):
#                     dest_q, dest_s = working_rc.right_root_at(dest_row, dest_col)
#                     if (dest_q, dest_s) == (a, b):
#                         # Place it here
#                         working_rc = working_rc.toggle_ref_at(dest_row, dest_col)
#                         break
#             else:
#                 continue
#             break
    
#     return working_rc
#     """
#     Inverse algorithm based on Kogan-Kumar paper.
    
#     Starting from top row, look right-to-left for reflections to remove.
#     After removing, look to the right for restoration configurations.
#     """
#     if len(reflections) == 0:
#         return rc, ()
    
#     # Build dictionaries tracking reflections
#     pair_dict = {}
#     pair_dict_rev = {}
#     for a, b in reflections:
#         pair_dict[a] = pair_dict.get(a, set())
#         pair_dict[a].add(b)
#         pair_dict_rev[b] = a
    
#     working_rc = rc
#     rows = []
    
#     # Process rows from top to bottom
#     for row in range(1, len(rc) + 1):
#         # Look from right to left in this row
#         col = working_rc.cols + 15
#         while col >= 1:
#             if working_rc.has_element(row, col):
#                 q, s = working_rc.right_root_at(row, col)
                
#                 # Check if this is a reflection to remove (Figure 9 configuration)
#                 removed = False
#                 target_a = None
#                 target_b = None
                
#                 # Case 1: Direct match (s in pair_dict_rev and q == pair_dict_rev[s])
#                 if s in pair_dict_rev and q == pair_dict_rev[s]:
#                     target_a = q
#                     target_b = s
#                     removed = True
                    
#                 # Case 2: Chained reflection (both q and s map to same parent)
#                 elif (s in pair_dict_rev and q in pair_dict_rev and 
#                       pair_dict_rev[s] == pair_dict_rev[q]):
#                     target_a = pair_dict_rev[s]
#                     target_b = s
#                     removed = True
                
#                 if removed:
#                     # Remove the intersection
#                     working_rc = working_rc.toggle_ref_at(row, col)
                    
#                     # Update dictionaries
#                     pair_dict[target_a].discard(target_b)
#                     if len(pair_dict[target_a]) == 0:
#                         del pair_dict[target_a]
#                     if target_b in pair_dict_rev:
#                         del pair_dict_rev[target_b]
                    
#                     rows.append(row)
                    
#                     # Look to the RIGHT for Figure 10 restoration configurations
#                     for col2 in range(col + 1, working_rc.cols + 20):
#                         if not working_rc.has_element(row, col2):
#                             q2, s2 = working_rc.right_root_at(row, col2)
                            
#                             if q2 > s2:
#                                 continue
                            
#                             # Figure 10 case 1: (q, s) where q <= d < s and s not already tracked
#                             if q2 <= d < s2 and s2 not in pair_dict_rev:
#                                 working_rc = working_rc.toggle_ref_at(row, col2)
#                                 pair_dict_rev[s2] = q2
#                                 pair_dict[q2] = pair_dict.get(q2, set())
#                                 pair_dict[q2].add(s2)
#                                 break
                                
#                             # Figure 10 case 2: (ai, q) where s2 already tracked
#                             elif q2 > d and s2 in pair_dict_rev and q2 not in pair_dict_rev:
#                                 working_rc = working_rc.toggle_ref_at(row, col2)
#                                 parent = pair_dict_rev[s2]
#                                 pair_dict_rev[q2] = parent
#                                 pair_dict[parent].add(q2)
#                                 break
                    
#                     # Continue looking left from the start of the row
#                     col = working_rc.cols + 15
#                 else:
#                     col -= 1
#             else:
#                 col -= 1
    
#     if len(pair_dict) > 0:
#         raise AssertionError(f"Not all reflections removed: {pair_dict}")
    
#     return working_rc, tuple(reversed(rows))

# def reverse_kk_insert_from_reflections(rc: RCGraph, d, reflections):
#     """
#     Reverse Kogan-Kumar insertion using reflections.
    
#     Key insight: forward insertion can push elements DOWN to lower rows during rectification.
#     Reverse must pull elements back UP from lower rows.
#     """
#     if len(reflections) == 0:
#         return rc, ()
    
#     pair_dict = {}
#     pair_dict_rev = {}
#     for a, b in reflections:
#         pair_dict[a] = pair_dict.get(a, set())
#         pair_dict[a].add(b)
#         pair_dict_rev[b] = a
    
#     working_rc = rc
#     rows = []
    
#     # Process all rows looking for reflections to remove
#     for row in range(1, len(rc) + 1):
#         col = 1
#         while col <= working_rc.cols + 10:
#             if working_rc.has_element(row, col):
#                 q, s = working_rc.right_root_at(row, col)
                
#                 removed = False
#                 # Case 1: Direct match (s in pair_dict_rev and q == pair_dict_rev[s])
#                 if s in pair_dict_rev and q == pair_dict_rev[s]:
#                     working_rc = working_rc.toggle_ref_at(row, col)
#                     del pair_dict_rev[s]
#                     pair_dict[q].discard(s)
#                     if len(pair_dict[q]) == 0:
#                         del pair_dict[q]
#                     removed = True
#                     rows.append(row)
                
#                 # Case 2: Chained reflection
#                 elif s in pair_dict_rev and q in pair_dict_rev and pair_dict_rev[s] == pair_dict_rev[q]:
#                     working_rc = working_rc.toggle_ref_at(row, col)
#                     parent = pair_dict_rev[s]
#                     pair_dict[parent].discard(s)
#                     if len(pair_dict[parent]) == 0:
#                         del pair_dict[parent]
#                     del pair_dict_rev[s]
#                     removed = True
#                     rows.append(row)
                
#                 # After removing, restore bumped elements
#                 if removed:
#                     restored = False
                    
#                     # First, try looking to the right in same row
#                     for col2 in range(col + 1, working_rc.cols + 15):
#                         if not working_rc.has_element(row, col2):
#                             q2, s2 = working_rc.right_root_at(row, col2)
#                             if q2 > s2:
#                                 continue
                            
#                             # Restore row-root reflections
#                             if q2 <= d < s2 and s2 not in pair_dict_rev:
#                                 working_rc = working_rc.toggle_ref_at(row, col2)
#                                 pair_dict_rev[s2] = q2
#                                 pair_dict[q2] = pair_dict.get(q2, set())
#                                 pair_dict[q2].add(s2)
#                                 restored = True
#                                 break
#                             # Restore chained reflections
#                             elif d < q2 and s2 in pair_dict_rev and q2 not in pair_dict_rev:
#                                 working_rc = working_rc.toggle_ref_at(row, col2)
#                                 pair_dict_rev[q2] = pair_dict_rev[s2]
#                                 pair_dict[pair_dict_rev[s2]].add(q2)
#                                 restored = True
#                                 break
                    
#                     # If nothing to the right, look in rows below
#                     if not restored:
#                         for check_row in range(row + 1, len(working_rc) + 1):
#                             for check_col in range(1, working_rc.cols + 15):
#                                 if working_rc.has_element(check_row, check_col):
#                                     q_below, s_below = working_rc.right_root_at(check_row, check_col)
                                    
#                                     # Move this element back up if it's not in our reflection list
#                                     if (q_below, s_below) not in [(a, b) for a in pair_dict for b in pair_dict[a]]:
#                                         if s_below not in pair_dict_rev and q_below not in pair_dict_rev:
#                                             # This element doesn't belong to inserted reflections
#                                             # Try to find where to put it in the current row
#                                             for target_col in range(1, working_rc.cols + 15):
#                                                 if not working_rc.has_element(row, target_col):
#                                                     q_target, s_target = working_rc.right_root_at(row, target_col)
                                                    
#                                                     # Check if (q_below, s_below) matches (q_target, s_target)
#                                                     if q_below == q_target and s_below == s_target:
#                                                         # Move it!
#                                                         working_rc = working_rc.toggle_ref_at(check_row, check_col)
#                                                         working_rc = working_rc.toggle_ref_at(row, target_col)
#                                                         restored = True
#                                                         break
                                                
#                                                 if restored:
#                                                     break
                                        
#                                         if restored:
#                                             break
                                
#                                 if restored:
#                                     break
                            
#                             if restored:
#                                 break
                    
#                     # Start over from beginning of this row
#                     col = 1
#                 else:
#                     col += 1
#             else:
#                 col += 1
    
#     if len(pair_dict) > 0:
#         raise AssertionError(f"Not all reflections used up: {pair_dict} on {working_rc}")
    
#     return working_rc, tuple(reversed(rows))

# Test

def anti_crystal(rc):
    rc = rc.resize(len(rc.perm) - 1)
    rc_t_list = [list(row) for row in rc.transpose(len(rc))]
    result_list = [[] * (len(rc))]
    for col_num in range(len(rc)):
        try:
            result_list[col_num].append(len(rc.perm) - rc_t_list[col_num].pop())
        except IndexError:
            pass
    return RCGraph([tuple(row) for row in result_list]).normalize()

def seq_product(seq1, seq2):
    result = RCGraph()
    n = len(seq1)
    for i in range(n):
        result = result.resize(i+1).squash_product(seq1[i].resize(i+1)).squash_product(seq2[i].resize(i+1))
    return result

def rc_from_seq(seq):
    result_rc = RCGraph([()])
    for i, rc in enumerate(seq):
        result_rc = result_rc.resize(i+1).squash_product(rc.resize(i+1)).resize(len(seq))
    return result_rc

# def grass_tableau_rep(grass_rc):

def factorize(rc_graph):
    rc_graph = rc_graph.normalize()
    if len(rc_graph) == 1:
        return [rc_graph]
    if len(rc_graph.perm.descents()) <= 1:
        return [rc_graph]

if __name__ == "__main__":
    #from schubmult.utils.schub_lib import elem_sym_perms_q
    import sys
    import itertools

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)

    grass_perms = {}
    for p in perms:
        if len(p.descents()) <= 1 and max(p.trimcode, default=0) == 1:
            grass_perms[len(p.trimcode)] = grass_perms.get(len(p.trimcode), [])
            grass_perms[len(p.trimcode)].append(p)
    # for i in range(1, n):
    #     rc_map2 = {rc.resize(i): seq + (RCGraph().resize(i),) for rc, seq in rc_map.items()}
    #     for gp in grass_perms[i]:
            
    #         for rc in RCGraph.all_rc_graphs(gp, i):
    for perm in perms:
        if perm.inv == 0:
            continue
        for rc in RCGraph.all_rc_graphs(perm, len(perm.trimcode)):
            # for k in  range(len(perm.trimcode), n):
            for perm2 in perms:
                if perm2.inv == 0:
                    continue
                if max(perm2.descents(), default=0) <= min(perm.descents(), default=0):
                    bottom_size = min(perm.descents()) + 1
                    for rc2 in RCGraph.all_rc_graphs(perm2, bottom_size):
                        rc_mul = RCGraph([*rc2.squash_product(rc.vertical_cut(bottom_size)[0]), *rc[bottom_size:]])
                        # tab = RootTableau.from_rc_graph(rc_g).weight_tableau
                        # print(f"Testing {perm}, {gp}, {rc=}, {rc_g=} {tab=}")
                        # rc_mul = rc.resize(k)
                        # for letter in reversed(tab.row_word):
                        #     rc_mul = rc_mul.squash_product(RCGraph(tuple([()]*(letter-1)) + ((k,),)).resize(k))
                        assert (Sx(perm) * Sx(perm2)).get(rc_mul.perm, 0) != S.Zero, f"Dammit on {perm}, {perm2}, {rc=}, {rc2=} {rc_mul=}"

                
                
    #     rc_map = rc_map2
    # print(rc_map)
    # print(len(rc_map))
    # print(sum(len(RCGraph.all_rc_graphs(perm, len(perm.trimcode))) for perm in perms))
    # for rc1, rc2 in itertools.product(rc_map.keys(), repeat=2):
    #     seq1 = rc_map[rc1]
    #     seq2 = rc_map[rc2]
    #     perm1 = rc1.perm
    #     perm2 = rc2.perm
    #     if len(perm1.trimcode) < n-1 or len(perm2.trimcode) < n-1:
    #         continue
    #     prod_rc = seq_product(seq1, seq2)
    #     #prod_rc = rc_from_seq(prod_seq)
    #     assert (Sx(perm1)*Sx(perm2)).get(prod_rc.perm, 0) != S.Zero, f"Failed on {perm1}, {perm2}, {rc1=} {rc2=} {prod_rc=} {seq1=} {seq2=}"
    #     print("Yay pgieons")

    