from schubmult import Permutation, elem_sym_perms, elem_sym_perms_q, q_vector, phi_d, tau_d, RCGraph

def reverse_kk_insert(rc: RCGraph, d, reflections):
    pair_dict = {}
    pair_dict_rev = {}
    for a,b in reflections:
        pair_dict[a] = pair_dict.get(a, [])
        pair_dict[a].append(b)
        pair_dict_rev[b] = a
    
    working_rc = rc
    rows = []
    for row in range(1, len(rc) + 1):
        for col in range(working_rc.cols + 5, 0, -1):
            if working_rc.has_element(row, col):
                q, s = working_rc.right_root_at(row, col)
                look_right = False
                if s in pair_dict_rev and q == pair_dict_rev[s]:
                    working_rc = working_rc.toggle_ref_at(row, col)
                    del pair_dict_rev[s]
                    pair_dict[q].remove(s)
                    if len(pair_dict[q]) == 0:
                        del pair_dict[q]
                elif s in pair_dict_rev and q in pair_dict_rev and pair_dict_rev[s] == pair_dict_rev[q]:
                    working_rc = working_rc.toggle_ref_at(row, col)
                    pair_dict[pair_dict_rev[s]].remove(s)
                    if len(pair_dict[pair_dict_rev[s]]) == 0:
                        del pair_dict[pair_dict_rev[s]]
                    del pair_dict_rev[s]

                    rows.append(row)
                else:
                    continue
                found = False
                for col2 in range(col + 1, working_rc.cols + 1):
                    if not working_rc.has_element(row, col2):
                        q2, s2 = working_rc.right_root_at(row, col2)
                        if q2 > s2:
                            continue
                        if q2 <= d and s2 > d and s2 not in pair_dict_rev:
                            working_rc = working_rc.toggle_ref_at(row, col2)
                            pair_dict_rev[q2] = s2
                            pair_dict[q2] = pair_dict.get(q2, [])
                            pair_dict[q2].append(s2)
                            found = True
                            break
                        elif d < q2 and s2 in pair_dict_rev and q2 not in pair_dict_rev:
                            working_rc = working_rc.toggle_ref_at(row, col2)
                            pair_dict_rev[q2] = pair_dict_rev[s2]
                            pair_dict[pair_dict_rev[s2]].append(q2)
                            found = True
                            break
                if not found:
                    rows.append(row)
    assert len(pair_dict) == 0, f"Not all reflections used up: {pair_dict} {working_rc}"
    return working_rc, tuple(reversed(rows))

# Test
if __name__ == "__main__":
    #from schubmult.utils.schub_lib import elem_sym_perms_q
    import sys

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)

    for perm in perms:
        rc: RCGraph
        for rc in RCGraph.all_rc_graphs(perm, n - 1):
            for d in range(1, n):
                for r in range(1, d + 1):
                    rc2, reflections = rc.kogan_kumar_insert(d, [r], return_reflections=True)
                    print(f"Testing perm {perm}, d={d}, r={r}, rc={rc}, rc2={rc2}, reflections={reflections}")
                    rc_reversed, rows = reverse_kk_insert(rc2, d, reflections)
                    assert rc_reversed == rc, f"Failed on perm {perm}, d={d}, r={r}, rc={rc} rc2={rc2} rc_reversed={rc_reversed}"
                    print("Success")

    