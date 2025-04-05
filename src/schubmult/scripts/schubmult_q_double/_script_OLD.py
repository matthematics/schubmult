import sys
from functools import cached_property

import numpy as np
from symengine import expand, symarray, sympify

from schubmult._base_argparse import schub_argparse
from schubmult.perm_lib import (
    check_blocks,
    code,
    count_less_than,
    inv,
    inverse,
    is_parabolic,
    longest_element,
    medium_theta,
    mulperm,
    omega,
    permtrim,
    q_var,
    q_vector,
    reduce_q_coeff,
    trimcode,
    uncode,
)
from schubmult.schubmult_double import compute_positive_rep, div_diff, posify
from schubmult.schubmult_q_double._funcs import (
    factor_out_q_keep_factored,
    # mult_poly,
    nil_hecke,
    schubmult,
    schubmult_db,
)


class _gvars:
    @cached_property
    def n(self):
        return 100

    # @cached_property
    # def fvar(self):
    #     return 100

    @cached_property
    def var1(self):
        return tuple(symarray("x", self.n).tolist())

    @cached_property
    def var2(self):
        return tuple(symarray("y", self.n).tolist())

    @cached_property
    def var3(self):
        return tuple(symarray("z", self.n).tolist())

    @cached_property
    def var_r(self):
        return symarray("r", 100)


_vars = _gvars()


def _display_full(coeff_dict, args, formatter, posified=None, var2=_vars.var2, var3=_vars.var3):
    raw_result_dict = {}
    mult = args.mult

    perms = args.perms

    ascode = args.ascode
    same = args.same
    check = args.check
    msg = args.msg
    display_positive = args.display_positive
    expa = args.expa
    slow = args.slow
    nilhecke_apply = False
    subs_dict2 = {}
    for i in range(1, 100):
        sm = var2[1]
        for j in range(1, i):
            sm += _vars.var_r[j]
        subs_dict2[var2[i]] = sm

    coeff_perms = list(coeff_dict.keys())
    coeff_perms.sort(key=lambda x: (inv(x), *x))

    for perm in coeff_perms:
        val = coeff_dict[perm]
        if expand(val) != 0:
            try:
                int(val)
            except Exception:
                val2 = 0
                if display_positive and not posified:
                    q_dict = factor_out_q_keep_factored(val)
                    for q_part in q_dict:
                        try:
                            val2 += q_part * int(q_dict[q_part])
                        except Exception:
                            if same:
                                to_add = q_part * expand(sympify(q_dict[q_part]).xreplace(subs_dict2))
                                val2 += to_add
                            else:
                                try:
                                    if len(perms) == 2:
                                        u = tuple(permtrim([*perms[0]]))
                                        v = tuple(permtrim([*perms[1]]))
                                    if (
                                        len(perms) == 2
                                        and code(inverse(perms[1])) == medium_theta(inverse(perms[1]))
                                        and not mult
                                        and not slow
                                        and not nilhecke_apply
                                    ):
                                        val2 += q_part * q_dict[q_part]
                                    else:
                                        q_part2 = q_part
                                        if not mult and not nilhecke_apply and len(perms) == 2:
                                            qv = q_vector(q_part)
                                            u2, v2, w2 = u, v, perm
                                            u2, v2, w2, qv, did_one = reduce_q_coeff(u2, v2, w2, qv)
                                            while did_one:
                                                u2, v2, w2, qv, did_one = reduce_q_coeff(u2, v2, w2, qv)
                                            q_part2 = np.prod(
                                                [q_var[i + 1] ** qv[i] for i in range(len(qv))],
                                            )
                                            if q_part2 == 1:
                                                # reduced to classical coefficient
                                                val2 += q_part * posify(
                                                    q_dict[q_part],
                                                    u2,
                                                    v2,
                                                    w2,
                                                    var2,
                                                    var3,
                                                    msg,
                                                    False,
                                                )
                                            else:
                                                val2 += q_part * compute_positive_rep(
                                                    q_dict[q_part],
                                                    var2,
                                                    var3,
                                                    msg,
                                                    False,
                                                )
                                        else:
                                            val2 += q_part * compute_positive_rep(
                                                q_dict[q_part],
                                                var2,
                                                var3,
                                                msg,
                                                False,
                                            )
                                except Exception as e:
                                    if mult:
                                        print(
                                            "warning; --display-positive is on but result is not positive",
                                            file=sys.stderr,
                                        )
                                        val2 = val
                                        break
                                    print(
                                        f"error; write to schubmult@gmail.com with the case {perms=} {perm=} {val=} {coeff_dict.get(perm,0)=}",
                                    )
                                    print(f"Exception: {e}")
                                    import traceback

                                    traceback.print_exc()
                                    exit(1)
                    if not same and check and expand(val - val2) != 0:
                        if mult:
                            val2 = val
                        else:
                            print(
                                f"error: value not equal; write to schubmult@gmail.com with the case {perms=} {perm=} {val2=} {coeff_dict.get(perm,0)=}",
                            )
                            exit(1)
                    val = val2
            if expa:
                val = expand(val)
            if val != 0:
                if ascode:
                    raw_result_dict[tuple(trimcode(perm))] = val
                    if formatter:
                        print(f"{trimcode(perm)!s}  {formatter(val)}")
                else:
                    raw_result_dict[tuple(perm)] = val
                    if formatter:
                        print(f"{perm!s}  {formatter(val)}")
    return raw_result_dict


def main(argv=None):
    if argv is None:
        argv = sys.argv
    var2 = tuple(symarray("y", 100))
    var3 = tuple(symarray("z", 100))
    try:
        sys.setrecursionlimit(1000000)

        args, formatter = schub_argparse(
            "schubmult_q_double",
            "Compute coefficients of products of quantum double Schubert polynomials in the same or different sets of coefficient variables",
            yz=True,
            quantum=True,
            argv=argv[1:],
        )
        subs_dict2 = {}
        for i in range(1, 100):
            sm = var2[1]
            for j in range(1, i):
                sm += _vars.var_r[j]
            subs_dict2[var2[i]] = sm

        mult = args.mult  # noqa: F841
        mulstring = args.mulstring  # noqa: F841

        perms = args.perms

        ascode = args.ascode
        msg = args.msg
        display_positive = args.display_positive
        pr = args.pr
        parabolic_index = [int(s) for s in args.parabolic]
        parabolic = len(parabolic_index) != 0
        slow = args.slow
        nil_N = 0
        nilhecke = False
        nilhecke_apply = False
        same = args.same
        if same:
            var3 = var2

        if args.nilhecke is not None:
            nilhecke = True
            nil_N = args.nilhecke
        if args.nilhecke_apply is not None:
            nil_N = args.nilhecke_apply
            nilhecke_apply = True

        if ascode:
            for i in range(len(perms)):
                perms[i] = tuple(permtrim(uncode(perms[i])))
        else:
            for i in range(len(perms)):
                if len(perms[i]) < 2 and (len(perms[i]) == 0 or perms[i][0] == 1):
                    perms[i] = (1, 2)
                perms[i] = tuple(permtrim([*perms[i]]))

        if nilhecke:
            coeff_dict = nil_hecke({(1, 2): 1}, perms[0], nil_N)
        elif nilhecke_apply:
            coeff_dict0 = nil_hecke({(1, 2): 1}, perms[0], nil_N, var2, var2)
            coeff_dict = {(1, 2): 0}
            for v in coeff_dict0:
                coeff_dict[(1, 2)] += coeff_dict0[v] * div_diff(v, perms[1], var2, var3)
        else:
            coeff_dict = {perms[0]: 1}
            for perm in perms[1:]:
                if not slow:
                    coeff_dict = schubmult_db(coeff_dict, perm, var2, var3)
                else:
                    coeff_dict = schubmult(coeff_dict, perm, var2, var3)
                # if mult:
                #     for v in var2:
                #         globals()[str(v)] = v
                #     for v in var3:
                #         globals()[str(v)] = v
                #     for v in var_x:
                #         globals()[str(v)] = v
                #     for v in q_var:
                #         globals()[str(v)] = v

                # mul_exp = eval(mulstring)
                # coeff_dict = mult_poly(coeff_dict, mul_exp)

        posified = False
        if parabolic:
            if display_positive:
                posified = True
            w_P = longest_element(parabolic_index)
            w_P_prime = [1, 2]
            coeff_dict_update = {}
            for w_1 in coeff_dict:
                val = coeff_dict[w_1]
                q_dict = factor_out_q_keep_factored(val)
                for q_part in q_dict:
                    qv = q_vector(q_part)
                    w = [*w_1]
                    good = True
                    parabolic_index2 = []
                    for i in range(len(parabolic_index)):
                        if omega(parabolic_index[i], qv) == 0:
                            parabolic_index2 += [parabolic_index[i]]
                        elif omega(parabolic_index[i], qv) != -1:
                            good = False
                            break
                    if not good:
                        continue
                    w_P_prime = longest_element(parabolic_index2)
                    if not check_blocks(qv, parabolic_index):
                        continue
                    w = permtrim(mulperm(mulperm(w, w_P_prime), w_P))
                    if not is_parabolic(w, parabolic_index):
                        continue

                    w = tuple(permtrim(w))

                    new_q_part = np.prod(
                        [
                            q_var[index + 1 - count_less_than(parabolic_index, index + 1)] ** qv[index]
                            for index in range(len(qv))
                            if index + 1 not in parabolic_index
                        ],
                    )

                    try:
                        new_q_part = int(new_q_part)
                    except Exception:
                        pass
                    q_val_part = q_dict[q_part]
                    if display_positive:
                        try:
                            q_val_part = int(q_val_part)
                        except Exception:
                            if same:
                                q_val_part = expand(sympify(q_val_part).xreplace(subs_dict2))
                            else:
                                try:
                                    if len(perms) == 2 and q_part == 1:
                                        u = permtrim([*perms[0]])
                                        v = permtrim([*perms[1]])
                                        q_val_part = posify(
                                            q_dict[q_part],
                                            tuple(u),
                                            tuple(v),
                                            w_1,
                                            var2,
                                            var3,
                                            msg,
                                            False,
                                        )
                                    else:
                                        qv = q_vector(q_part)
                                        u2, v2, w2 = perms[0], perms[1], w_1
                                        u2, v2, w2, qv, did_one = reduce_q_coeff(u2, v2, w2, qv)
                                        while did_one:
                                            u2, v2, w2, qv, did_one = reduce_q_coeff(u2, v2, w2, qv)
                                        q_part2 = np.prod(
                                            [q_var[i + 1] ** qv[i] for i in range(len(qv))],
                                        )
                                        if q_part2 == 1:
                                            q_val_part = posify(
                                                q_dict[q_part],
                                                u2,
                                                v2,
                                                w2,
                                                var2,
                                                var3,
                                                msg,
                                                False,
                                            )
                                        else:
                                            q_val_part = compute_positive_rep(
                                                q_dict[q_part],
                                                var2,
                                                var3,
                                                msg,
                                                False,
                                            )
                                except Exception as e:
                                    print(
                                        f"error; write to schubmult@gmail.com with the case {perms=} {perm=} {q_part*q_val_part=} {coeff_dict.get(w_1,0)=}",
                                    )
                                    print(f"Exception: {e}")
                                    exit(1)
                            coeff_dict_update[w] = coeff_dict_update.get(w, 0) + new_q_part * q_val_part

            coeff_dict = coeff_dict_update

        raw_result_dict = {}
        if pr or formatter is None:
            raw_result_dict = _display_full(coeff_dict, args, formatter, posified)
        if formatter is None:
            return raw_result_dict
    except BrokenPipeError:
        pass


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv))
