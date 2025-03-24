import pytest
from ast import literal_eval
import re


from schubmult._tests import get_json


def check_positive(v2, same, subs_dict2):
    # if same, should be no minus signs
    from symengine import sympify, expand
    from schubmult.schubmult_double import compute_positive_rep
    from schubmult.schubmult_double._funcs import _vars
    from schubmult.schubmult_q_double import factor_out_q_keep_factored

    var2 = _vars.var2
    var3 = _vars.var3
    import sys
    q_dict = factor_out_q_keep_factored(v2)
    for k, v in q_dict.items():
        if same:
            v = expand(sympify(v).xreplace(subs_dict2))
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
                assert expand(sym_val) == expand(
                    compute_positive_rep(sym_val, var2, var3, do_pos_neg=False)
                )
    return True


def assert_dict_good(
    v_tuple,
    input_dict,
    ret_dict,
    same=True,
    display_positive=False,
    slow = False
):
    # print(f"{input_dict=}")
    
    from schubmult.schubmult_q_double import schubmult_db, schubmult
    from symengine import symarray, expand, sympify

    var_a = symarray("y", 100).tolist()
    var_b = symarray("z", 100).tolist()
    var_q = symarray("q", 100).tolist()
    var_r = symarray("r", 100).tolist()
    #print(f"{ret_dict=}",file=sys.stderr)
    subs_dict2 = {}
    for i in range(1, 100):
        sm = var_a[1]
        for j in range(1, i):
            sm += var_r[j]
        subs_dict2[var_a[i]] = sm
    if slow:
        coeff_dict = schubmult(input_dict, v_tuple, var2=var_a, var3=var_a if same else var_b,q_var=var_q)
    else:
        coeff_dict = schubmult_db(input_dict, v_tuple, var2=var_a, var3=var_a if same else var_b,q_var=var_q)

    #print(f"{coeff_dict=}",file=sys.stderr)
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
        coeff_dict = {k: expand(sympify(v).xreplace(subs_dict2)) for k, v in coeff_dict.items()}
    for k, v in coeff_dict.items():
        import sys        
        if k in ret_dict:
            if display_positive and not same:
                print(f"boofle {same=} {display_positive=} {v=} {ret_dict[k]=}", file=sys.stderr)
            assert (expand(v) == expand(ret_dict[k]))
        else:
            assert expand(v) == 0
    for k in ret_dict.keys():
        if display_positive:
            assert check_positive(ret_dict[k], same, subs_dict2)
        v = coeff_dict[k]
        # v = sympify(v).xreplace(subs_dict2)
        assert expand(v) == expand(ret_dict[k])


def parse_ret(lines, ascode, unformat):
    from schubmult.perm_lib import uncode

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
        ret_dict[(tuple(literal_eval(k)) if not ascode else tuple(uncode(literal_eval(k))))] = v
    return ret_dict



base_dir = "script_tests/data/schubmult_q_double"

json_files_data_args = [
    "test_gen1",
    "test_mixed_positive",    
    "test_same_positive",    
]


@pytest.mark.parametrize("json_file", json_files_data_args)
def test_with_same_args_exec(capsys, json_file):
    from schubmult.perm_lib import uncode, permtrim

    args = get_json(f"{base_dir}/{json_file}")
    # print(f"{json_file=} {args=} input_data")
    from schubmult.schubmult_q_double._script import main

    mult = args["mult"]  # noqa: F841
    mulstring = args["mulstring"]  # noqa: F841

    perms = args["perms"]

    ascode = args["ascode"]
    same = args["same"]
    msg = args["msg"]  # noqa: F841
    display_positive = args["display_positive"]  # noqa: F841
    pr = args["pr"]  # noqa: F841
    disp_mode = args["disp_mode"]
    slow = args["slow"]

    from latex2sympy2_extended import latex2sympy
    from symengine import sympify

    unformat = {"basic": lambda v: sympify(v), "latex": lambda v: sympify(latex2sympy(v))}
    # print(f"{args=}")
    # print(f"{args['cmd_line']=}")
    main(args["cmd_line"])
    lines = capsys.readouterr()
    lines = str(lines.out).split("\n")

    ret_dict = parse_ret(lines, ascode, unformat[disp_mode])
    v_tuple = (tuple(perms[1]) if not ascode else tuple(uncode(perms[1])))
    input_dict = {tuple(permtrim(perms[0])) if not ascode else tuple(permtrim(uncode(perms[0]))): 1}
    
    assert_dict_good(v_tuple, input_dict, ret_dict, same, display_positive, slow)


if __name__ == "__main__":
    # test_coeff_equal_exec()
    test_with_same_args_exec()
