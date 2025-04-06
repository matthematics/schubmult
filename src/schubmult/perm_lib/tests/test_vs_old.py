def test():
    pass

# from schubmult.schub_lib.tests.legacy_perm_lib import (
#     code,
#     count_bruhat,
#     count_less_than,
#     cycle,
#     dominates,
#     get_cycles,
#     has_bruhat_ascent,
#     has_bruhat_descent,
#     inv,
#     inverse,
#     medium_theta,
#     mu_A,
#     mulperm,
#     omega,
#     one_dominates,
#     p_trans,
#     permtrim,
#     phi1,
#     sg,
#     strict_theta,
#     theta,
#     trimcode,
#     uncode,
#     elem_sym_func_q,
#     compute_vpathdicts,
#     double_elem_sym_q
# )

# from schubmult.schub_lib.quantum_double import _vars
# import schubmult.perm_lib.perm_lib as pl
# import schubmult.schub_lib.schub_lib as sl
# import schubmult.poly_lib.poly_lib as sply

# def test_schubmult_db():
#     perm_dict = {(3,4,1,2): 1}
#     perm_dict_N = {pl.Permutation(k): v for k,v in perm_dict.items()}
#     v = (4,3,1,2)
#     v_N = pl.Permutation(v)
#     var2=_vars.var2
#     var3=_vars.var3
#     q_var=_vars.q_var

#     th = medium_theta(inverse(v))
#     th_N = pl.medium_theta(~v_N)
#     assert list(th) == list(th_N)
#     while th[-1] == 0:
#         th.pop()
#     while th_N[-1] == 0:
#         th_N.pop()
#     assert th == th_N
#     mu = permtrim(uncode(th))
#     mu_N = pl.uncode(th_N)
#     assert list(mu) == list(mu_N)
#     vmu = permtrim(mulperm([*v], mu))
#     vmu_N = v_N*mu_N
#     assert list(vmu) == list(vmu_N)
#     inv_vmu = inv(vmu)
#     inv_vmu_N = vmu_N.inv
#     assert inv_vmu == inv_vmu_N
#     inv_mu = inv(mu)
#     inv_mu_N = mu_N.inv
#     assert inv_mu == inv_mu_N

#     ret_dict = {}
#     ret_dict_N = {}
#     thL = len(th)
#     vpathdicts = compute_vpathdicts(th, vmu, True)
#     vpathdicts_N = sl.compute_vpathdicts(th_N, vmu_N, True)
#     for u, val in perm_dict.items():
#         u_N = pl.Permutation(u)
#         inv_u = inv(u)
#         vpathsums = {u: {(1, 2): val}}
#         vpathsums_N = {u_N: {pl.Permutation((1, 2)): val}}
#         for index in range(thL):
#             if index > 0 and th[index - 1] == th[index]:
#                 continue
#             mx_th = 0
#             for vp in vpathdicts[index]:
#                 for v2, vdiff, s in vpathdicts[index][vp]:
#                     mx_th = max(mx_th, th[index] - vdiff)
#             mx_th_N = 0
#             for vp_N in vpathdicts_N[index]:
#                 for v2_N, vdiff_N, s_N in vpathdicts_N[index][vp_N]:
#                     mx_th_N = max(mx_th_N, th_N[index] - vdiff_N)
#             assert mx_th_N == mx_th
#             if index < len(th) - 1 and th[index] == th[index + 1]:
#                 mx_th1 = 0
#                 for vp in vpathdicts[index + 1]:
#                     for v2, vdiff, s in vpathdicts[index + 1][vp]:
#                         mx_th1 = max(mx_th1, th[index + 1] - vdiff)
#                 newpathsums = {}
#                 newpathsums_N = {}
#                 for up in vpathsums:
#                     up_N = pl.Permutation(up)
#                     newpathsums0 = {}
#                     newpathsums0_N = {}
#                     inv_up = inv(up)
#                     newperms = double_elem_sym_q(up, mx_th, mx_th1, th[index], q_var)
#                     newperms_N = sl.double_elem_sym_q(up_N, mx_th, mx_th1, th[index], q_var)
#                     for vv in vpathdicts[index]:
#                         sumval = vpathsums[up].get(vv, 0)
#                         vv_N = pl.Permutation(vv)
#                         sumval_N = vpathsums_N[up_N].get(vv_N, 0)
#                         assert sumval == sumval_N
#                         if sumval == 0:
#                             continue
#                         for v2, vdiff2, s2 in vpathdicts[index][vv]:
#                             v2_N = pl.Permutation(v2)
#                             assert (v2_N, vdiff2, s2) in vpathdicts_N[index][vv_N]
#                             for up1, udiff1, mul_val1 in newperms:
#                                 up1_N = pl.Permutation(up1)
#                                 assert (up1_N, udiff1, mul_val1) in newperms_N
#                                 esim1 = (
#                                     elem_sym_func_q(
#                                         th[index],
#                                         index + 1,
#                                         up,
#                                         up1,
#                                         vv,
#                                         v2,
#                                         udiff1,
#                                         vdiff2,
#                                         var2,
#                                         var3,
#                                     )
#                                     * mul_val1
#                                     * s2
#                                 )
#                                 esim1_N = (
#                                     sply.elem_sym_func_q(
#                                         th[index],
#                                         index + 1,
#                                         up_N,
#                                         up1_N,
#                                         vv_N,
#                                         v2_N,
#                                         udiff1,
#                                         vdiff2,
#                                         var2,
#                                         var3,
#                                     )
#                                     * mul_val1
#                                     * s2
#                                 )
#                                 assert esim1_N == esim1
#                                 mulfac = sumval * esim1
#                                 if (up1, udiff1, mul_val1) not in newpathsums0:
#                                     assert (up1_N, udiff1, mul_val1) not in newpathsums0_N
#                                     newpathsums0[(up1, udiff1, mul_val1)] = {}
#                                     newpathsums0_N[(up1_N, udiff1, mul_val1)] = {}
#                                 # newpathsums0[(up1, udiff1, mul_val1
#                                 newpathsums0[(up1, udiff1, mul_val1)][v2] = newpathsums0[(up1, udiff1, mul_val1)].get(v2, 0) + mulfac
#                                 newpathsums0[(up1_N, udiff1, mul_val1)][v2_N] = newpathsums0_N[(up1_N, udiff1, mul_val1)].get(v2_N, 0) + mulfac
#             else:
#                 newpathsums = {}
#                 for up in vpathsums:
#                     inv_up = inv(up)
#                     newperms = elem_sym_perms_q(
#                         up,
#                         min(mx_th, (inv_mu - (inv_up - inv_u)) - inv_vmu),
#                         th[index],
#                         q_var,
#                     )
#                     for up2, udiff, mul_val in newperms:
#                         if up2 not in newpathsums:
#                             newpathsums[up2] = {}
#                         for v in vpathdicts[index]:
#                             sumval = vpathsums[up].get(v, 0) * mul_val
#                             if sumval == 0:
#                                 continue
#                             for v2, vdiff, s in vpathdicts[index][v]:
#                                 newpathsums[up2][v2] = newpathsums[up2].get(
#                                     v2,
#                                     0,
#                                 ) + s * sumval * elem_sym_func_q(
#                                     th[index],
#                                     index + 1,
#                                     up,
#                                     up2,
#                                     v,
#                                     v2,
#                                     udiff,
#                                     vdiff,
#                                     var2,
#                                     var3,
#                                 )
#             vpathsums = newpathsums
#         toget = tuple(vmu)
#         ret_dict = add_perm_dict({ep: vpathsums[ep].get(toget, 0) for ep in vpathsums}, ret_dict)
    #return ret_dict
