import re
from ast import literal_eval

import pytest

from schubmult.utils.test_utils import get_json, load_json_test_names
from schubmult.utils.parsing import parse_coeff
from schubmult import GeneratingSet

def check_positive(v, coprod, same, var_r):
    # if same, should be no minus signs
    from schubmult.symbolic import expand, sympify

    from schubmult import compute_positive_rep
    var2 = GeneratingSet("y")
    var3 = GeneratingSet("z")
    import sys
    # print(f"{var2=} {var3=}",file=sys.stderr)

    if same:
        if coprod:
            subs_dict = {var_r[i]: -var_r[i] for i in range(1, 100)}
            v = expand(v.subs(subs_dict))
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

    #print(f"{vals=}", file=sys.stderr)
    for val in vals:
        if val.find("(") != -1:
            val += ")"
        try:
            if int(val) >= 0:
                continue
        except ValueError:
            pass
        sym_val = expand(parse_coeff(val))
        #print(f"{sym_val=} not an int", file=sys.stderr)
        if coprod:
            assert expand(sym_val) == expand(
                compute_positive_rep(sym_val, var3, var2),
            )
        else:
            assert expand(sym_val) == expand(
                compute_positive_rep(sym_val, var2, var3),
            )
    #print(f"Donebaby {coprod=}", file=sys.stderr)
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
    from schubmult import Permutation
    from schubmult.symbolic import expand, sympify
    import sys
    from schubmult import schub_coprod_double, schubmult_double

    var_a = GeneratingSet("y")
    var_b = GeneratingSet("z")

    if coprod:
        coeff_dict = schub_coprod_double(v_tuple, indices, var2=var_a, var3=var_a if same else var_b)
    else:
        coeff_dict = schubmult_double(input_dict, v_tuple, var2=var_a, var3=var_a if same else var_b)
    # for k in ret_dict.keys():
    #     # print(f"{k=} {type(k)=} {type(k[0])=} {type(k[1])=} {type(ret_dict[k])=}")
    #print(f"{coeff_dict=}",file=sys.stderr)
    # print(f"{ret_dict=}")
    var_r = GeneratingSet("r")
    subs_dict2 = {}
    for i in range(1, 100):
        sm = var_a[1]
        for j in range(1, i):
            sm += var_r[j]
        subs_dict2[var_a[i]] = sm
    # if same and display_positive:
    #     coeff_dict = {k: expand(sympify(v).subs(subs_dict2)) for k, v in coeff_dict.items()}
    for k, v in coeff_dict.items():
        if expand(v) == 0:
            assert (k not in ret_dict) or (expand(v) == expand(ret_dict[k]))
        else:
            #print(f"{k=} {type(k[0])=} {type(k[1])=} {coeff_dict[k]=}")
            assert expand(v) == expand(ret_dict[k])
    for k in ret_dict.keys():
        if display_positive and not same:
            assert check_positive(ret_dict[k], coprod, same, var_r)
        v = coeff_dict[k]
        # v = sympify(v).xreplace(subs_dict2)
        assert expand(v) == expand(ret_dict[k])


def parse_ret(lines, ascode, coprod, unformat):
    import sys

    from schubmult import permtrim, uncode
    from schubmult import Permutation

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
            ret_dict[(Permutation(literal_eval(k)) if not ascode else (uncode(literal_eval(k))))] = v
    else:
        for line in lines:
            line = str(line)
            charo = "[" if ascode else "("
            charc = "]" if ascode else ")"
            first_split = "\\) +\\(" if not ascode else "\\] +\\["
            second_split = ")  " if not ascode else "]  "
            jn = ")," if not ascode else "],"
            try:
                s, vf = re.split(first_split, line, maxsplit=1)
                f, v = vf.split(second_split, maxsplit=1)  # re.split(second_split, vf)
                evlaf = f"({charo}{f + jn + s}{charc})"
                k1, k2 = literal_eval(evlaf)
                if ascode:
                    k1 = uncode(k1)
                    k2 = uncode(k2)
                    # print("uncoded")
                # print(f"{k1=} {k2=}")
                k = (Permutation(k2), Permutation(k1))
            except Exception as e:
                # print(f"boingfish {line=} {e=}", file=sys.stderr)
                continue
            try:
                v = unformat(v)                
            except Exception as e:
                # print(f"bingfish {line=} {v=} {e=}", file=sys.stderr)
                continue            
            ret_dict[k] = v
            #print(f"ret_dict[{k=}] = {v} {type(k)=} {type(k[0])=} {type(ret_dict[k])=}")
    return ret_dict


#
base_dir = "schubmult_double"

json_files_data_args = load_json_test_names(base_dir)


@pytest.mark.parametrize("json_file", json_files_data_args)
def test_with_same_args_exec(capsys, json_file):
    from schubmult import permtrim, uncode

    args = get_json(f"{base_dir}/{json_file}")
    # print(f"{json_file=} {args=} input_data")
    from schubmult._scripts.schubmult_double import main
    from schubmult import Permutation

    mult = args["mult"]  # noqa: F841
    mulstring = args["mulstring"]  # noqa: F841

    perms = args["perms"]

    ascode = args["ascode"]
    coprod = args["coprod"]
    same = args["same"]
    msg = args["msg"]  # noqa: F841
    down = args["down"]  # noqa: F841
    display_positive = args["display_positive"]
    pr = args["pr"]  # noqa: F841
    disp_mode = args["disp_mode"]

    from latex2sympy2_extended import latex2sympy
    from schubmult.symbolic import sympify

    unformat = {"basic": lambda v: parse_coeff(v), "latex": lambda v: parse_coeff(str(latex2sympy(v)))}
    # print(f"{args=}")
    # print(f"{args['cmd_line']=}")
    ret_dict = main(args["cmd_line"])
    lines = capsys.readouterr()
    lines = str(lines.out).split("\n")

    if disp_mode != "raw":
        ret_dict = parse_ret(lines, ascode, coprod, unformat[disp_mode])
    elif coprod:
        if ascode:
            ret_dict = {((uncode(list(k[0]))), (uncode(list(k[1])))): v for k, v in ret_dict.items()}
        else:
            ret_dict = {((Permutation(k[0])), (Permutation(list(k[1])))): v for k, v in ret_dict.items()}
    elif ascode:
        ret_dict = {uncode(list(k)): v for k, v in ret_dict.items()}
    v_tuple = (Permutation(perms[1]) if not ascode else uncode(perms[1])) if not coprod else (Permutation(perms[0]) if not ascode else (uncode(perms[0])))
    input_dict = {(permtrim(perms[0])) if not ascode else (permtrim(uncode(perms[0]))): 1}
    indices = tuple(perms[1])    
    assert_dict_good(v_tuple, input_dict, ret_dict, coprod, indices, same, display_positive)


if __name__ == "__main__":
    # test_coeff_equal_exec()
    class fope:
        pass
    spoingle = fope()
    st = ""
    with(open("floil","r")) as f:
        st = str(f.read())
    boing = fope()
    boing.out = st
    spoingle.readouterr = lambda: boing
    test_with_same_args_exec(spoingle,
                            "schubmult_double_mixed-var_display-positive_coprod_code_1_0_2_0_1T2_4_display-mode_latex")
