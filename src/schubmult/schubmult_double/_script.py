import sys
from functools import cached_property

import sympy
from symengine import expand, sympify

from schubmult._base_argparse import schub_argparse
from schubmult.logging import get_logger
from schubmult.perm_lib import (
    Permutation,
    add_perm_dict,
    code,
    inv,
    mu_A,
    permtrim,
    split_perms,
    theta,
    trimcode,
    uncode,
)
from schubmult.poly_lib import GeneratingSet, efficient_subs, is_indexed
from schubmult.schub_lib import (
    will_formula_work,
)

# from schubmult.schubmult_double._vars import var_x, var, var_r
from schubmult.schubmult_double._funcs import (
    # !TEMP
    # compute_positive_rep,
    mult_poly,
    mult_poly_down,
    posify,
    schub_coprod,
    schubmult,
    schubmult_down,
)

logger = get_logger(__name__)


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


_vars = _gvars()


def _display(val):
    print(val)


subs_dict = {}
for i in range(1, 100):
    sm = _vars.var2[1]
    for j in range(1, i):
        sm += _vars.var_r[j]
    subs_dict[_vars.var2[i]] = sm


def pre_posify(perms, perm, val, check, check_val, same, down, var2, var3, msg, subs_dict):
    try:
        return int(val)
    except Exception:
        if same:
            val = efficient_subs(sympify(val), subs_dict).expand()  # expand(sympify(val).xreplace(subs_dict))
        else:
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
            # except Exception:
            #     _display(
            #         f"error; write to schubmult@gmail.com with the case {perms=} {perm=} {val=} {check_val=}",
            #     )
            #     exit(1)
            # if check and expand(val - check_coeff_dict.get(perm, 0)) != 0:
            if check and expand(val - check_val) != 0:
                _display(
                    f"error; write to schubmult@gmail.com with the case {perms=} {perm=} {val=} {check_val=}",
                )
                logger.debug("Yep it's here")
                exit(1)
    return val

def flip_symbol_signs(val):
    subs_dict={}
    for s in val.free_symbols:
        if is_indexed(s):
            subs_dict[s] = -s
    return val.subs(subs_dict)


def _display_full(
    coeff_dict,
    args,
    formatter,
    var2,
    var3,
    check=True,
):
    subs_dict2 = {}
    for i in range(1, 100):
        sm = var2[1]
        for j in range(1, i):
            sm += _vars.var_r[j]
        subs_dict2[var2[i]] = sm
    raw_result_dict = {}
    perms = args.perms
    ascode = args.ascode
    coprod = args.coprod
    msg = args.msg
    same = args.same
    display_positive = args.display_positive
    perms[0] = Permutation(perms[0])
    pos = list(perms[1])
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

    mu_W = uncode(theta(~perms[0]))

    the_top_perm = perms[0]*mu_W

    muA = uncode(mu_A(code(mu_W), poso))
    muB = uncode(mu_A(code(mu_W), pos2))

    coeff_perms = list(coeff_dict.keys())
    if coprod:
        perm_pairs = coeff_perms
        if ascode:
            width = max(
                [len(str(trimcode(perm[0])) + " " + str(trimcode(perm[1]))) for perm in perm_pairs],
            )
        else:
            width = max([len(str(perm[0]) + " " + str(perm[1])) for perm in perm_pairs])

        for firstperm, secondperm in perm_pairs:
            val = coeff_dict[(firstperm,secondperm)]
            if same and display_positive:
                val = efficient_subs(sympify(val), subs_dict2).expand()
            if val != 0:
                if display_positive and not same:
                    if val != 0:
                        try:
                            val = int(expand(val))
                        except Exception as e:
                            val2 = posify(
                                flip_symbol_signs(val),
                                firstperm * muA,
                                secondperm * muB,
                                the_top_perm,
                                var2,
                                var3,
                                msg,
                                False,
                            )
                            val2 = flip_symbol_signs(val2)
                            if check and expand(val - val2) != 0:
                                _display(
                                    f"error; write to schubmult@gmail.com with the case {perms=}\n{code(firstperm)=} {code(secondperm)=}\n{val2=}\n{val=}",
                                )
                                _display(f"{firstperm*muA=} {secondperm*muB=} {the_top_perm=}")
                                exit(1)
                            val = val2
                    else:
                        val = 0
                if val != 0:
                    if not ascode:
                        width2 = width - len(str(permtrim(firstperm))) - len(str(permtrim(secondperm)))
                        raw_result_dict[(permtrim(firstperm), Permutation(secondperm))] = val
                        if formatter:
                            _display(
                                f"{permtrim(firstperm)}{' ':>{width2}}{Permutation(secondperm)}  {formatter(val)}",
                            )
                    else:
                        width2 = width - len(str(trimcode(firstperm))) - len(str(trimcode(secondperm)))
                        raw_result_dict[(tuple(trimcode(firstperm)), tuple(trimcode(secondperm)))] = val
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
            # if val != 0:
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
        var2 = GeneratingSet("y")
        var3 = GeneratingSet("z")
        sys.setrecursionlimit(1000000)

        # TEMP
        sympy.init_printing()

        args, formatter = schub_argparse(
            "schubmult_double",
            "Compute coefficients of products of double Schubert polynomials in the same or different sets of coefficient variables",
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
        check = args.check
        display_positive = args.display_positive
        pr = args.pr
        # debug = args.debug

        # logger.log(logging.DEBUG, f"main boing 1 {var2=}{var3=}{same=}")
        if same:
            var3 = var2
        # logger.log(logging.DEBUG, f"main boing 2 {var2=}{var3=}{same=}")
        posified = False
        if coprod:
            if ascode:
                perms[0] = uncode(perms[0])
            pos = [*perms[1]]
            pos.sort()
            mperm = Permutation(perms[0])

            coeff_dict = schub_coprod(mperm, pos, var2, var3)

            if pr or formatter is None:
                # logger.log(logging.DEBUG, f"main {var2=}{var3=}{same=}")
                raw_result_dict = _display_full(
                    coeff_dict,
                    args,
                    formatter,
                    var2=var2,
                    var3=var3,
                )
            if formatter is None:
                return raw_result_dict
        else:
            if ascode:
                for i in range(len(perms)):
                    perms[i] = uncode(perms[i])
            else:
                for i in range(len(perms)):
                    if len(perms[i]) < 2 and (len(perms[i]) == 0 or perms[i][0] == 1):
                        perms[i] = (1, 2)
                    perms[i] = Permutation(perms[i])

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
            # preprocess positivity
            if display_positive and len(perms) == 2 and will_formula_work(perms[0], perms[1]) and not mult and not down and not same:
                coeff_dict = {}
                th = theta(perms[1])
                muv = uncode(th)
                muvn1v = (~muv) * perms[1]
                coeff_dict2 = {perms[0]: 1}
                coeff_dict2 = schubmult(coeff_dict2, muv, var2, var3)
                for perm, val in coeff_dict2.items():
                    w = perm * muvn1v
                    if inv(w) + inv(muvn1v) == inv(perm):
                        coeff_dict[Permutation(w)] = val
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
                                coeff_dict4[w],
                                u,
                                perm,
                                w,
                                var2,
                                var3,
                                msg,
                            )
                        coeff_dict3 = add_perm_dict(coeff_dict4, coeff_dict3)
                    coeff_dict2 = coeff_dict3
                coeff_dict = coeff_dict2
                posified = True
            elif not posified:
                coeff_dict = check_coeff_dict

            if not posified and display_positive:
                # print(f"{coeff_dict=}")
                coeff_dict = {k: pre_posify(perms, k, v, check, check_coeff_dict.get(k, 0), same, down, var2, var3, msg, subs_dict) for k, v in coeff_dict.items()}

            if pr or formatter is None:
                raw_result_dict = _display_full(
                    coeff_dict,
                    args,
                    formatter,
                    var2,
                    var3,
                )

            if formatter is None:
                return raw_result_dict
    except BrokenPipeError:
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
