import schubmult.poly_lib.variables as schub_vars
import schubmult.schub_lib.schub_lib as schub_lib
import schubmult.perm_lib as pl
import schubmult.schub_lib.tests.legacy_perm_lib as opl
q_var = schub_vars.GeneratingSet("q")



def test_double_elem_sym_q():
    ret_list = {}
    ret_list_N = {}
    u= tuple([7,5,1,6,4,2,3])
    u_N = pl.Permutation(u)
    p1 = 3
    p2 = 5
    k = 6

    perms1 = opl.elem_sym_perms_q(u, p1, k, q_var)
    perms1_N = schub_lib.elem_sym_perms_q(u_N, p1, k, q_var)
    iu = opl.inverse(u)
    iu_N = ~u_N
    for perm1, udiff1, mul_val1 in perms1:
        perm1_N = pl.Permutation(perm1)
        assert (perm1_N, udiff1, mul_val1) in perms1_N, f"Permutation {perm1} not found in new result"
        perms2 = opl.elem_sym_perms_q(perm1, p2, k, q_var)
        perms2_N = schub_lib.elem_sym_perms_q(perm1_N, p2, k, q_var)
        cycles1 = opl.get_cycles(opl.mulperm([*iu],[*perm1]))
        cycles1_N = (iu_N*perm1_N).get_cycles()
        assert cycles1 == cycles1_N
        cycles1_dict = {}
        for c in cycles1:
            if c[-1] not in cycles1_dict:
                cycles1_dict[c[-1]] = []
            cycles1_dict[c[-1]] += [set(c)]
        ip1 = opl.inverse(perm1)
        ip1_N = ~perm1_N
        assert ip1 == list(ip1_N)
        for perm2, udiff2, mul_val2 in perms2:
            perm2_N = pl.Permutation(perm2)
            assert (perm2_N, udiff2, mul_val2) in perms2_N, f"Permutation {perm2} not found in new result"
            cycles2 = opl.get_cycles(opl.mulperm([*ip1],[*perm2]))
            cycles2_N = (ip1_N*perm2_N).get_cycles()
            assert cycles2 == cycles2_N
            good = True
            for i in range(len(cycles2)):
                c2 = cycles2[i]
                if c2[-1] not in cycles1_dict:
                    continue
                for c1_s in cycles1_dict[c2[-1]]:
                    for a in range(len(c2) - 2, -1, -1):
                        if c2[a] in c1_s:
                            good = False
                            break
                    if not good:
                        break
                if not good:
                    break

            if good:
                if (perm1, udiff1, mul_val1) not in ret_list:
                    ret_list[(perm1, udiff1, mul_val1)] = []
                ret_list[(perm1, udiff1, mul_val1)] += [(perm2, udiff2, mul_val2)]
    return ret_list



def test_elem_sym_q():
    # Test the function with a sample permutation
    old_perm = tuple([5,1,4,2,3])
    orig_perm = pl.Permutation(old_perm)
    p = 3
    k = 5
    result_old = opl.elem_sym_perms_q(old_perm, p, k)
    result_new = schub_lib.elem_sym_perms_q(orig_perm, p, k)
    assert len(result_old) == len(result_new), "Length of results do not match"
    for i in range(len(result_old)):
        # print(f"{result_old[i]=} {result_new[i]=}")
        assert pl.Permutation(result_old[i][0]) == result_new[i][0], f"Permutation mismatch at index {i}"
        assert result_old[i][2] == result_new[i][2], f"Value mismatch at index {i}"
    # def double_elem_sym_q(u, p1, p2, k, q_var=q_var):
    old_perm = tuple([7,5,1,6,4,2,3])
    orig_perm = pl.Permutation(old_perm)
    p1 = 3
    p2 = 5
    k = 6
    result_old = opl.double_elem_sym_q(old_perm, p1, p2, k)
    result_new = schub_lib.double_elem_sym_q(orig_perm, p1, p2, k)
    assert len(result_old) == len(result_new), "Length of results do not match"
    print("OLD")
    print("\n".join([str((str(b),str(result_old[b]))) for b in result_old]))
    print("NEW")
    print("\n".join([str((str(b),str(result_new[b]))) for b in result_new]))
    for k, v in result_old.items():
        k2 = (pl.Permutation(k[0]),k[1],k[2])
        assert k2 in result_new, f"Key {k2} not found in new result"
        v2 = result_new[k2]
        v = sorted(v)
        v2 = sorted(v2)
        assert len(v) == len(v2)
        for i in range(len(v)):
            assert pl.Permutation(v[i][0]) == v2[i][0]
            assert v[i][2] == v2[i][2]
        
        #double_elem_sym_q
        # assert pl.Permutation(result_old[i][1]) == result_new[i][1], f"Value mismatch at index {i}"
        # assert pl.Permutation(result_old[i][2]) == result_new[i][2], f"Last j mismatch at index {i}"
    # Print the result
    # for perm, val, last_j in result:
    #     # print(f"Permutation: {perm}, Value: {val}, Last j: {last_j}")