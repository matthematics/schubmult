import sys
from functools import cached_property

import numpy as np
import sympy
from symengine import expand, symarray, sympify

from schubmult._base_argparse import schub_argparse
from schubmult.perm_lib import (
    add_perm_dict,
    code,
    inv,
    inverse,
    mu_A,
    mulperm,
    permtrim,
    theta,
    trimcode,
    uncode,
    will_formula_work,
)

# from schubmult.schubmult_double._vars import var_x, var, var_r
from schubmult.schubmult_double._funcs import (
    compute_positive_rep,
    mult_poly,
    mult_poly_down,
    posify,
    schubmult,
    schubmult_down,
    split_perms,
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


def _display(val):
    print(val)


def _display_full(
    coeff_dict,
    args,
    formatter,
    var2,
    var3,
    posified=None,
    check_coeff_dict=None,
    kperm=None,
    N=None,
):
    subs_dict2 = {}
    for i in range(1, 100):
        sm = var2[1]
        for j in range(1, i):
            sm += _vars.var_r[j]
        subs_dict2[var2[i]] = sm
    raw_result_dict = {}
    perms = args.perms
    mult = args.mult
    ascode = args.ascode
    coprod = args.coprod
    check = args.check
    msg = args.msg
    down = args.down
    same = args.same
    display_positive = args.display_positive

    coeff_perms = list(coeff_dict.keys())
    if coprod:
        pos = perms[1]
        pos2 = []
        last_descent = -1
        poso = []
        for i in range(len(perms[0]) - 1):
            if perms[0][i] > perms[0][i + 1]:
                last_descent = i + 1
        for i in range(1, last_descent + 1):
            if i not in pos:
                pos2 += [i - 1]
            else:
                poso += [i - 1]

        mu_W = uncode(theta(inverse(perms[0])))

        the_top_perm = tuple(permtrim(mulperm(list(perms[0]), mu_W)))

        muA = uncode(mu_A(code(mu_W), poso))
        muB = uncode(mu_A(code(mu_W), pos2))
        subs_dict = {}
        inv_kperm = inv(kperm)
        inverse_kperm = inverse(kperm)
        var2neg = np.array([-var2[i] for i in range(len(var2))])
        var3neg = np.array([-var3[i] for i in range(len(var3))])

        for i in range(1, 100):
            if i <= N:
                subs_dict[_vars.var1[i]] = var2[i]
            else:
                subs_dict[_vars.var1[i]] = var3[i - N]

        coeff_perms.sort(key=lambda x: (inv(x), *x))

        perm_pairs = []

        for perm in coeff_perms:
            downperm = mulperm(list(perm), inverse_kperm)
            if inv(downperm) == inv(perm) - inv_kperm:
                flag = True
                for i in range(N):
                    if downperm[i] > N:
                        flag = False
                        break
                if not flag:
                    continue
                firstperm = downperm[0:N]
                secondperm = [downperm[i] - N for i in range(N, len(downperm))]
                perm_pairs += [[permtrim(firstperm), permtrim(secondperm)]]

        if ascode:
            width = max(
                [len(str(trimcode(perm[0])) + " " + str(trimcode(perm[1]))) for perm in perm_pairs],
            )
        else:
            width = max([len(str(perm[0]) + " " + str(perm[1])) for perm in perm_pairs])


        for perm in coeff_perms:
            val = coeff_dict[perm]
            downperm = mulperm(list(perm), inverse_kperm)
            if inv(downperm) == inv(perm) - inv_kperm:
                flag = True
                for i in range(N):
                    if downperm[i] > N:
                        flag = False
                        break
                if not flag:
                    continue
                firstperm = downperm[0:N]
                secondperm = [downperm[i] - N for i in range(N, len(downperm))]
                val = sympify(val).subs(subs_dict)

                if same and display_positive:
                    val = expand(sympify(val).subs(subs_dict2))

                if val != 0:
                    if display_positive and not same:
                        if val != 0:
                            val2 = posify(
                                val,
                                tuple(permtrim(mulperm(firstperm, muA))),
                                tuple(permtrim(mulperm(secondperm, muB))),
                                the_top_perm,
                                tuple(var2neg.tolist()),
                                tuple(var3neg.tolist()),
                                msg,
                                False,
                            )
                            if expand(val - val2) != 0:
                                _display(
                                    f"error; write to schubmult@gmail.com with the case {perms=}\n{code(firstperm)=} {code(secondperm)=}\n{val2=}\n{val=}",
                                )
                                _display(
                                    f"{code(tuple(permtrim(mulperm(firstperm,muA))))=},{code(tuple(permtrim(mulperm(secondperm,muB))))=},{code(the_top_perm)=}\n{expand(val-val2)=}",
                                )
                                exit(1)
                            val = val2
                        else:
                            val = 0
                    if val != 0:
                        if not ascode:
                            width2 = (
                                width
                                - len(str(permtrim(firstperm)))
                                - len(str(permtrim(secondperm)))
                            )
                            raw_result_dict[
                                (tuple(permtrim(firstperm)), tuple(permtrim(secondperm)))
                            ] = val
                            if formatter:
                                _display(
                                    f"{tuple(permtrim(firstperm))}{' ':>{width2}}{tuple(permtrim(secondperm))}  {formatter(val)}",
                                )
                        else:
                            width2 = (
                                width
                                - len(str(trimcode(firstperm)))
                                - len(str(trimcode(secondperm)))
                            )
                            raw_result_dict[
                                (tuple(trimcode(firstperm)), tuple(trimcode(secondperm)))
                            ] = val
                            if formatter:
                                _display(
                                    f"{trimcode(firstperm)}{' ':>{width2}}{trimcode(secondperm)}  {formatter(val)}",
                                )
    else:
        if ascode:
            width = max([len(str(trimcode(perm))) for perm in coeff_dict.keys()])
        else:
            width = max([len(str(perm)) for perm in coeff_dict.keys()])

        coeff_perms = list(coeff_dict.keys())
        coeff_perms.sort(key=lambda x: (inv(x), *x))

        for perm in coeff_perms:
            val = coeff_dict[perm]
            if val != 0:
                notint = False
                try:
                    int(val)
                except Exception:
                    notint = True
                if notint and display_positive:
                    if same:
                        subs_dict = {}
                        for i in range(1, 100):
                            sm = var2[1]
                            for j in range(1, i):
                                sm += _vars.var_r[j]
                            subs_dict[var2[i]] = sm
                        val = expand(sympify(coeff_dict[perm]).xreplace(subs_dict))
                    else:
                        try:
                            if len(perms) == 2 and not posified and not mult:
                                if not down:
                                    val = posify(
                                        val,
                                        perms[0],
                                        perms[1],
                                        perm,
                                        var2,
                                        var3,
                                        msg,
                                    )
                                else:
                                    val = posify(
                                        val,
                                        perm,
                                        perms[1],
                                        perms[0],
                                        var2,
                                        var3,
                                        msg,
                                    )
                            elif not posified and not mult:
                                val = compute_positive_rep(val, var2, var3, msg)
                            elif not posified:
                                val = compute_positive_rep(val, var2, var3, msg, do_pos_neg=False)
                        except Exception:
                            if mult:
                                _display(
                                    "warning; --display-positive is on but result is not positive",
                                    file=sys.stderr,
                                )
                            else:
                                _display(
                                    f"error; write to schubmult@gmail.com with the case {perms=} {perm=} {val=} {check_coeff_dict.get(perm,0)=}",
                                )
                                exit(1)
                        if check and expand(val - check_coeff_dict.get(perm, 0)) != 0:
                            _display(
                                f"error; write to schubmult@gmail.com with the case {perms=} {perm=} {val=} {check_coeff_dict.get(perm,0)=}",
                            )
                            exit(1)
                if val != 0:
                    if ascode:
                        raw_result_dict[tuple(trimcode(perm))] = val
                        if formatter:
                            _display(f"{trimcode(perm)!s:>{width}}  {formatter(val)}")
                    else:
                        raw_result_dict[tuple(perm)] = val
                        if formatter:
                            _display(f"{perm!s:>{width}}  {formatter(val)}")
    return raw_result_dict


def main(argv=None):
    if argv is None:
        argv = sys.argv

    try:
        var2 = tuple(symarray("y", 100).tolist())
        var3 = tuple(symarray("z", 100).tolist())
        sys.setrecursionlimit(1000000)

        # TEMP
        sympy.init_printing()

        args, formatter = schub_argparse(
            "schubmult_double",
            "Compute coefficients of product of double Schubert polynomials in the same or different sets of coefficient variables",
            argv=argv[1:],
            yz=True,
        )

        mult = args.mult
        mulstring = args.mulstring

        perms = args.perms

        ascode = args.ascode
        coprod = args.coprod
        same = args.same
        msg = args.msg
        down = args.down
        display_positive = args.display_positive
        pr = args.pr

        # logger.log(logging.DEBUG, f"main boing 1 {var2=}{var3=}{same=}")
        if same:
            # logger.log(logging.DEBUG, f"main OOO {same=}")
            var3 = var2
        # logger.log(logging.DEBUG, f"main boing 2 {var2=}{var3=}{same=}")
        posified = False
        if coprod:
            if ascode:
                mperm = tuple(permtrim(uncode(perms[0])))
            else:
                mperm = tuple(permtrim(perms[0]))

            perms[0] = mperm
            pos = perms[1]

            k = len(pos)
            n = len(perms[0])
            kcd = [pos[i] - i - 1 for i in range(len(pos))] + [n + 1 - k for i in range(k, n)]
            N = len(kcd)

            kperm = inverse(uncode(kcd))
            coeff_dict = {tuple(kperm): 1}
            coeff_dict = schubmult(coeff_dict, perms[0], _vars.var1, var2)

            if pr or formatter is None:
                # logger.log(logging.DEBUG, f"main {var2=}{var3=}{same=}")
                raw_result_dict = _display_full(
                    coeff_dict,
                    args,
                    formatter,
                    posified=posified,
                    kperm=kperm,
                    var2=var2,
                    var3=var3,
                    N=N,
                )
            if formatter is None:
                return raw_result_dict
        else:
            if ascode:
                for i in range(len(perms)):
                    perms[i] = tuple(permtrim(uncode(perms[i])))
            else:
                for i in range(len(perms)):
                    if len(perms[i]) < 2 and (len(perms[i]) == 0 or perms[i][0] == 1):
                        perms[i] = (1, 2)
                    perms[i] = tuple(permtrim([*perms[i]]))

            size = 0
            orig_perms = [*perms]
            while len(perms) != size:
                size = len(perms)
                perms = split_perms(perms)

            coeff_dict = {perms[0]: 1}
            check_coeff_dict = {perms[0]: 1}

            # if mult:
            #     for v in var2:
            #         ()[str(v)] = v
            #     for v in var3:
            #         globals()[str(v)] = v
            #     for v in _vars.var1:
            #         globals()[str(v)] = v

            if down:
                for perm in orig_perms[1:]:
                    check_coeff_dict = schubmult_down(check_coeff_dict, perm, var2, var3)
                if mult:
                    mul_exp = eval(mulstring)
                    check_coeff_dict = mult_poly_down(check_coeff_dict, mul_exp)
            else:
                for perm in orig_perms[1:]:
                    check_coeff_dict = schubmult(check_coeff_dict, perm, var2, var3)
                # coeff_dict = check_coeff_dict
                if mult:
                    mul_exp = eval(mulstring)
                    check_coeff_dict = mult_poly(check_coeff_dict, mul_exp)

            if (
                display_positive
                and len(perms) == 2
                and will_formula_work(perms[0], perms[1])
                and not mult
                and not down
                and not same
            ):
                coeff_dict = {}
                th = theta(perms[1])
                muv = uncode(th)
                muvn1v = mulperm(inverse(muv), perms[1])
                coeff_dict2 = {perms[0]: 1}
                coeff_dict2 = schubmult(coeff_dict2, muv, var2, var3)
                for perm, val in coeff_dict2.items():
                    w = mulperm([*perm], muvn1v)
                    if inv(w) + inv(muvn1v) == inv(perm):
                        coeff_dict[tuple(permtrim(w))] = val
                posified = True

            if display_positive and len(perms) > 2 and not mult and not same:
                coeff_dict2 = dict(coeff_dict)
                for perm in perms[1:]:
                    coeff_dict3 = {}
                    for u in coeff_dict2:
                        coeff_dict4 = {u: 1}
                        coeff_dict4 = schubmult(coeff_dict4, perm, var2, var3)
                        for w in coeff_dict4:
                            coeff_dict4[w] = coeff_dict2[u] * posify(
                                coeff_dict4[w], u, perm, w, var2, var3, msg,
                            )
                        coeff_dict3 = add_perm_dict(coeff_dict4, coeff_dict3)
                    coeff_dict2 = coeff_dict3
                coeff_dict = coeff_dict2
                posified = True
            elif not posified:
                coeff_dict = check_coeff_dict


            if pr or formatter is None:
                raw_result_dict = _display_full(
                    coeff_dict,
                    args,
                    formatter,
                    var2,
                    var3,
                    posified=posified,
                    check_coeff_dict=check_coeff_dict,
                )

            if formatter is None:
                return raw_result_dict
    except BrokenPipeError:
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
