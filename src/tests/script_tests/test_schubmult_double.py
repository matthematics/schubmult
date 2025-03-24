import pytest
from ast import literal_eval
import re


from schubmult._tests import get_json, load_json_test_names


def check_positive(v, coprod, same, var_r):
    # if same, should be no minus signs
    from symengine import sympify, expand
    from schubmult.schubmult_double import compute_positive_rep
    from schubmult.schubmult_double._funcs import _vars

    var2 = _vars.var2
    var3 = _vars.var3
    import sys

    if same:
        if coprod:
            subs_dict = {var_r[i]: -var_r[i] for i in range(1, 100)}
            v = expand(v.subs(subs_dict))
        return str(v).find("-") == -1
    else:
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

        print(f"{vals=}", file=sys.stderr)
        for val in vals:
            if val.find("(") != -1:
                val += ")"
            try:
                if int(val) >= 0:
                    continue
            except ValueError:
                pass
            sym_val = expand(sympify(val))
            print(f"{sym_val=} not an int", file=sys.stderr)
            if coprod:
                assert expand(sym_val) == expand(
                    compute_positive_rep(sym_val, var3, var2, do_pos_neg=False)
                )
            else:
                assert expand(sym_val) == expand(
                    compute_positive_rep(sym_val, var2, var3, do_pos_neg=False)
                )
    print(f"Donebaby {coprod=}", file=sys.stderr)
    return True


def assert_dict_good(
    v_tuple,
    input_dict,
    ret_dict,
    coprod,
    indices,
    same=True,
    display_positive=False,
):
    # print(f"{input_dict=}")
    from schubmult.schubmult_double import schubmult, schub_coprod
    from symengine import symarray, expand, sympify

    var_a = symarray("y", 100).tolist()
    var_b = symarray("z", 100).tolist()

    if coprod:
        coeff_dict = schub_coprod(v_tuple, indices, var2=var_a, var3=var_a if same else var_b)
    else:
        coeff_dict = schubmult(input_dict, v_tuple, var2=var_a, var3=var_a if same else var_b)

    # print(f"{coeff_dict=}")
    # print(f"{ret_dict=}")
    var_r = tuple(symarray("r", 100))
    subs_dict2 = {}
    for i in range(1, 100):
        sm = var_a[1]
        for j in range(1, i):
            sm += var_r[j]
        subs_dict2[var_a[i]] = sm
    if same and display_positive:
        coeff_dict = {k: expand(sympify(v).subs(subs_dict2)) for k, v in coeff_dict.items()}
    for k, v in coeff_dict.items():
        if expand(v) == 0:
            assert (k not in ret_dict) or (expand(v) == expand(ret_dict[k]))
        else:
            assert expand(v) == expand(ret_dict[k])
    for k in ret_dict.keys():
        if display_positive:
            assert check_positive(ret_dict[k], coprod, same, var_r)
        v = coeff_dict[k]
        # v = sympify(v).xreplace(subs_dict2)
        assert expand(v) == expand(ret_dict[k])


def parse_ret(lines, ascode, coprod, unformat):
    from schubmult.perm_lib import uncode, permtrim

    ret_dict = {}
    if not coprod:
        for line in lines:
            try:
                k, v = line.lstrip().split("  ")
            except Exception:
                continue
            try:
                v = unformat(v)
            except Exception:
                continue
            ret_dict[(tuple(literal_eval(k)) if not ascode else tuple(uncode(literal_eval(k))))] = v
    else:
        for line in lines:
            line = str(line)
            charo = "[" if ascode else "("
            charc = "]" if ascode else ")"
            first_split = "\\) +\\(" if not ascode else "\\] +\\["
            second_split = ")  " if not ascode else "]  "
            jn = ")," if not ascode else "],"
            try:
                s, vf = re.split(first_split, line)
                f, v = vf.split(second_split)  # re.split(second_split, vf)
                evlaf = f"({charo}{f + jn + s}{charc})"
                k1, k2 = literal_eval(evlaf)
                if ascode:
                    k1 = tuple(permtrim(uncode(k1)))
                    k2 = tuple(permtrim(uncode(k2)))
                k = (k2, k1)
            except Exception:
                continue
            try:
                v = unformat(v)
            except Exception:
                continue
            ret_dict[k] = v
    return ret_dict


# json_files_data = [
#     "schubmult_double_cf1",
#     "schubmult_double_cf2",
#     "schubmult_double_coprod_mixed",
#     "schubmult_double_coprod_mixed_positive",
# ]


# @pytest.mark.parametrize("json_file", json_files_data)
# def test_coeff_equal_exec(capsys, json_file):
#     import re
#     from schubmult.perm_lib import uncode, permtrim

#     input_data = get_json(json_file)
#     from schubmult.schubmult_double._script import main
#     from symengine import sympify, symarray, expand

#     # sys.argv = input_data["cmd_line"]
#     ascode = input_data.get("code", False)
#     # print(f"BOOB {json_file=} {input_data=}")
#     input_dict = {
#         literal_eval(k) if not ascode else tuple(permtrim(uncode(literal_eval(k)))): v
#         for k, v in input_data.get("input_dict",{}).items()
#     }
#     ret_dict = {}
#     coprod = input_data.get("coprod", False)
#     indices = input_data.get("indices", None)
#     same = not input_data.get("mixed", False)
#     main(input_data["cmd_line"])
#     lines = capsys.readouterr()
#     # print(f"BOOB {json_file=} {lines=} forple")
#     lines = str(lines.out).split("\n")

#     ret_dict = parse_ret(lines, ascode, coprod)
#     v_tuple = tuple(input_data["v_tuple"]) if not ascode else tuple(uncode(input_data["v_tuple"]))
#     # print("BOOB ASS")
#     assert_dict_good(v_tuple, input_dict, ret_dict, coprod, indices, same)
#     # print(f"BOOB {json_file=} yay")


base_dir = "schubmult_double"

json_files_data_args = load_json_test_names(base_dir)


@pytest.mark.parametrize("json_file", json_files_data_args)
def test_with_same_args_exec(capsys, json_file):
    from schubmult.perm_lib import uncode, permtrim

    args = get_json(f"{base_dir}/{json_file}")
    # print(f"{json_file=} {args=} input_data")
    from schubmult.schubmult_double._script import main

    mult = args["mult"]  # noqa: F841
    mulstring = args["mulstring"]  # noqa: F841

    perms = args["perms"]

    ascode = args["ascode"]
    coprod = args["coprod"]
    same = args["same"]
    msg = args["msg"]  # noqa: F841
    down = args["down"]  # noqa: F841
    display_positive = args["display_positive"]  # noqa: F841
    pr = args["pr"]  # noqa: F841
    disp_mode = args["disp_mode"]

    from latex2sympy2_extended import latex2sympy
    from symengine import sympify

    unformat = {"basic": lambda v: sympify(v), "latex": lambda v: sympify(latex2sympy(v))}
    # print(f"{args=}")
    # print(f"{args['cmd_line']=}")
    ret_dict = main(args["cmd_line"])
    lines = capsys.readouterr()
    lines = str(lines.out).split("\n")

    if disp_mode != "raw":
        ret_dict = parse_ret(lines, ascode, coprod, unformat[disp_mode])
    else:
        if ascode:
            ret_dict = {tuple(uncode(k)): v for k, v in ret_dict.items()}
    v_tuple = (
        (tuple(perms[1]) if not ascode else tuple(uncode(perms[1])))
        if not coprod
        else (tuple(perms[0]) if not ascode else tuple(uncode(perms[0])))
    )
    input_dict = {tuple(permtrim(perms[0])) if not ascode else tuple(permtrim(uncode(perms[0]))): 1}
    indices = tuple(perms[1])
    # print(f"{v_tuple=} {input_dict=} {indices=}")
    # print("BOOB ASS")
    assert_dict_good(v_tuple, input_dict, ret_dict, coprod, indices, same, display_positive)


if __name__ == "__main__":
    # test_coeff_equal_exec()
    test_with_same_args_exec()
