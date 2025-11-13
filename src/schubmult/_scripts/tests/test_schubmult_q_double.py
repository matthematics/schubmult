import re
from ast import literal_eval

import pytest

from schubmult.utils.test_utils import get_json, load_json_test_names
from schubmult.utils.parsing import parse_coeff
from schubmult import GeneratingSet, efficient_subs


def check_positive(v2, same, subs_dict2, var2, var3, q_var):
    # if same, should be no minus signs
    from schubmult.symbolic import expand, sympify

    from schubmult import compute_positive_rep
    
    from schubmult import factor_out_q_keep_factored
    

    q_dict = factor_out_q_keep_factored(v2)
    for k, v in q_dict.items():
        if same:
            v = expand(efficient_subs(parse_coeff(v),subs_dict2))
            return str(v).find("-") == -1
        try:
            if int(v) > 0:
                return True
        except Exception:
            pass
        work = str(v)
        work = re.sub(r"[)] [+] [(]|[)][*][(]", "|", work)
        work = re.sub(r"\s[(]|[)]", "", work)

        vals = [val for val in work.split("|")]

        for i in range(len(vals)):
            if vals[i].find("*(") == -1:
                vals[i] = re.sub(r"[)]|[(]", "", vals[i])

        # print(f"{vals=}", file=sys.stderr)
        for val in vals:
            if val.find("(") != -1:
                val += ")"
            try:
                if int(val) >= 0:
                    continue
            except ValueError:
                pass
            sym_val = expand(parse_coeff(val))
            # print(f"{sym_val=} not an int", file=sys.stderr)
            assert expand(sym_val) == expand(
                compute_positive_rep(sym_val, var2, var3),
            )
    return True


def assert_dict_good(v_tuple, input_dict, ret_dict, parabolic, same=True, display_positive=False, slow=False):
    # print(f"{input_dict=}")

    from schubmult.symbolic import expand, sympify

    from schubmult import schubmult_q_double, schubmult_q_double_fast, apply_peterson_woodward

    var_a = GeneratingSet("y")
    var_b = GeneratingSet("z")
    var_q = GeneratingSet("q")
    var_r = GeneratingSet("r")
    # print(f"{ret_dict=}",file=sys.stderr)
    subs_dict2 = {}
    for i in range(1):
        sm = var_a[1]
        for j in range(1, i):
            sm += var_r[j]
        subs_dict2[var_a[i]] = sm
    if slow:
        coeff_dict = schubmult_q_double(
            input_dict, v_tuple, var2=var_a, var3=var_a if same else var_b, q_var=var_q,
        )
    else:
        coeff_dict = schubmult_q_double_fast(
            input_dict, v_tuple, var2=var_a, var3=var_a if same else var_b, q_var=var_q,
        )
    if len(parabolic) > 0:
        parabolic_index = []
        start = 0
        # 1, 2 | 3
        for i in range(len(parabolic)):
            end = start + int(parabolic[i])
            parabolic_index += list(range(start + 1, end))
            # start += int(args.parabolic[i])
            start = end
        # [sum(int(args.parabolic[j]) for j in range(i+1)) for i in range(len(args.parabolic))]
        coeff_dict = apply_peterson_woodward(coeff_dict, parabolic_index)
    # print(f"{coeff_dict=}",file=sys.stderr)
    # print(f"{coeff_dict=}")
    # print(f"{ret_dict=}")
    var_r = GeneratingSet("r")
    subs_dict2 = {}
    for i in range(100):
        sm = var_a[1]
        for j in range(1, i):
            sm += var_r[j]
        subs_dict2[var_a[i]] = sm
    # if same and display_positive:
    #     coeff_dict = {k: expand(efficient_subs(sympify(v),subs_dict2)) for k, v in coeff_dict.items()}
    for k, v in coeff_dict.items():
        import sys

        if k in ret_dict:
            # if display_positive and not same:
            # print(f"boofle {same=} {display_positive=} {v=} {ret_dict[k]=}", file=sys.stderr)
            assert expand(v) == expand(ret_dict[k])
        else:
            # print(f"{k=} {expand(coeff_dict[k])=}", file=sys.stderr)
            assert expand(v) == 0
    for k in ret_dict.keys():
        if display_positive and not same:
            assert check_positive(ret_dict[k], same, subs_dict2, var_a, var_a if same else var_b, q_var=var_q)
        v = coeff_dict[k]
        # v = sympify(v).xreplace(subs_dict2)
        assert expand(v) == expand(ret_dict[k])


def parse_ret(lines, ascode, unformat):
    from schubmult import uncode, permtrim

    ret_dict = {}
    for line in lines:
        try:
            k, v = line.lstrip().split("  ")
        except Exception:
            continue
        try:
            v = unformat(v)
        except Exception:
            continue
        ret_dict[(permtrim(literal_eval(k)) if not ascode else (uncode(literal_eval(k))))] = v
    return ret_dict


base_dir = "schubmult_q_double"

json_files_data_args = load_json_test_names(base_dir)


@pytest.mark.parametrize("json_file", json_files_data_args)
def test_with_same_args_exec(capsys, json_file):
    from schubmult import permtrim, uncode


    args = get_json(f"{base_dir}/{json_file}")
    # print(f"{json_file=} {args=} input_data")
    from schubmult._scripts.schubmult_q_double import main

    mult = args["mult"]  # noqa: F841
    mulstring = args["mulstring"]  # noqa: F841

    perms = args["perms"]

    ascode = args["ascode"]
    same = args["same"]
    msg = args["msg"]  # noqa: F841
    display_positive = args["display_positive"]
    pr = args["pr"]  # noqa: F841
    disp_mode = args["disp_mode"]
    slow = args["slow"]
    parabolic = [int(bob) for bob in args['parabolic']]

    from latex2sympy2_extended import latex2sympy
    from schubmult.symbolic import sympify

    unformat = {"basic": lambda v: parse_coeff(v), "latex": lambda v: parse_coeff(str(latex2sympy(v)))}
    # print(f"{args=}")
    # print(f"{args['cmd_line']=}")
    ret_dict = main(args["cmd_line"])
    lines = capsys.readouterr()
    lines = str(lines.out).split("\n")

    if disp_mode != "raw":
        ret_dict = parse_ret(lines, ascode, unformat[disp_mode])
    elif ascode:
        ret_dict = {(uncode(k)): v for k, v in ret_dict.items()}
    v_tuple = permtrim(perms[1]) if not ascode else (uncode(perms[1]))
    input_dict = {(permtrim(perms[0])) if not ascode else (permtrim(uncode(perms[0]))): 1}

    assert_dict_good(v_tuple, input_dict, ret_dict, parabolic, same, display_positive, slow)


if __name__ == "__main__":
    # test_coeff_equal_exec()
    test_with_same_args_exec()
