import sys
from functools import cached_property

import numpy as np
import sympy
from symengine import sympify

from schubmult import GeneratingSet, div_diff, efficient_subs, q_vector
from schubmult.perm_lib import (
    Permutation,
    inv,
    longest_element,
    permtrim,
    uncode,
)
from schubmult.schub_lib.quantum_double import (
    factor_out_q_keep_factored,
    nil_hecke,
    q_posify,
    schubmult_q_double,
    schubmult_q_double_fast,
)
from schubmult.schub_lib.schub_lib import (
    check_blocks,
)
from schubmult.utils.argparse import schub_argparse
from schubmult.utils.perm_utils import (
    count_less_than,
    is_parabolic,
    omega,
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
        return GeneratingSet("x")

    @cached_property
    def var2(self):
        return GeneratingSet("y")

    @cached_property
    def var3(self):
        return GeneratingSet("z")

    @cached_property
    def var_r(self):
        return GeneratingSet("r")

    @cached_property
    def var_g1(self):
        return GeneratingSet("y")

    @cached_property
    def var_g2(self):
        return GeneratingSet("z")

    @cached_property
    def q_var(self):
        return GeneratingSet("q")


_vars = _gvars()

q_var = _vars.q_var

zero = sympify(0)

subs_dict2 = {}
for i in range(1, 100):
    sm = _vars.var2[1]
    for j in range(1, i):
        sm += _vars.var_r[j]
    subs_dict2[_vars.var2[i]] = sm


def _display_full(coeff_dict, args, formatter, var2=_vars.var2, var3=_vars.var3):  # noqa: ARG001
    ascode = args.ascode
    Permutation.print_as_code = ascode
    coeff_perms = list(coeff_dict.keys())
    coeff_perms.sort(key=lambda x: (inv(x), *x))

    raw_result_dict = {}

    for perm in coeff_perms:
        val = coeff_dict[perm]
        if val != 0:
            raw_result_dict[perm] = val
            if formatter:
                print(f"{sympy.sstr(perm)!s}  {formatter(val)}")
    return raw_result_dict


def sv_posify(val):
    # this has just y's, we want to rearrange
    # can we do this without an optimization
    val = sympify(sympy.simplify(efficient_subs(val, subs_dict2)))
    bingle_dict = {}
    for i in range(1, len(_vars.var_r) - 1):
        bingle_dict[_vars.var_r[i]] = _vars.var2[i + 1] - _vars.var2[i]  # sympy.Add(*[_vars.var2[i+1], - _vars.var2[i]],evaluate=False)
        # oh bay does that bar bangled banber bet bave space buckets of cheese
    # val = sympy.simplify(val)
    # return efficient_subs(val, bingle_dict)
    return val.xreplace(bingle_dict)


def main(argv=None):
    if argv is None:
        argv = sys.argv
    var2 = _vars.var2
    var3 = _vars.var3
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
        parabolic_index = []
        start = 0
        # 1, 2 | 3
        for i in range(len(args.parabolic)):
            end = start + int(args.parabolic[i])
            parabolic_index += list(range(start + 1, end))
            # start += int(args.parabolic[i])
            start = end
        # [sum(int(args.parabolic[j]) for j in range(i+1)) for i in range(len(args.parabolic))]
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
                perms[i] = uncode(perms[i])
        else:
            for i in range(len(perms)):
                if len(perms[i]) < 2 and (len(perms[i]) == 0 or perms[i][0] == 1):
                    perms[i] = Permutation([1, 2])
                perms[i] = permtrim(perms[i])

        if nilhecke:
            coeff_dict = nil_hecke({Permutation([1, 2]): 1}, perms[0], nil_N)
        elif nilhecke_apply:
            coeff_dict0 = nil_hecke({Permutation([1, 2]): 1}, perms[0], nil_N, var2, var2)
            coeff_dict = {Permutation([]): 0}
            for v in coeff_dict0:
                coeff_dict[Permutation([])] += coeff_dict0[v] * div_diff(v, perms[1], var2, var3)
        else:
            coeff_dict = {perms[0]: 1}
            for perm in perms[1:]:
                if not slow:
                    coeff_dict = schubmult_q_double_fast(coeff_dict, perm, var2, var3)
                else:
                    coeff_dict = schubmult_q_double(coeff_dict, perm, var2, var3)

        if parabolic:
            max_len = parabolic_index[-1] + 1
            # parabolic_index += list(range(parabolic_index[-1] + 2, max_len))
            w_P = longest_element(parabolic_index)
            # max_len = len(w_P)
            w_P_prime = Permutation([1, 2])
            coeff_dict_update = {}
            for w_1 in coeff_dict.keys():
                val = coeff_dict[w_1]
                q_dict = factor_out_q_keep_factored(val)
                for q_part in q_dict:
                    qv = q_vector(q_part)
                    w = w_1
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
                    w = (w * w_P_prime) * w_P
                    if not is_parabolic(w, parabolic_index):
                        continue

                    w = permtrim(w)
                    if len(w) > max_len:
                        continue
                    new_q_part = np.prod(
                        [q_var[index + 1 - count_less_than(parabolic_index, index + 1)] ** qv[index] for index in range(len(qv)) if index + 1 not in parabolic_index],
                    )
                    try:
                        new_q_part = int(new_q_part)
                    except Exception:
                        pass
                    q_val_part = q_dict[q_part]
                    coeff_dict_update[w] = coeff_dict_update.get(w, 0) + new_q_part * q_val_part
            coeff_dict = coeff_dict_update

        if display_positive and not nilhecke and not nilhecke_apply:
            if same:
                coeff_dict = {perm: sv_posify(val) for perm, val in coeff_dict.items()}
            else:
                coeff_dict = {perm: q_posify(perms[0], perms[1], perm, val, var2, var3, _vars.q_var, msg) for perm, val in coeff_dict.items()}

        raw_result_dict = {}

        if pr or formatter is None:
            raw_result_dict = _display_full(coeff_dict, args, formatter)
        if formatter is None:
            return raw_result_dict
    except BrokenPipeError:
        pass


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv))
