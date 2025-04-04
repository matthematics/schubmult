from ast import literal_eval

import pytest

from schubmult.utils import get_json, load_json_test_names


def assert_dict_good(v_tuple, input_dict, ret_dict):
    from schubmult.schub_lib.quantum import schubmult_q

    coeff_dict = schubmult_q(input_dict, v_tuple)
    for k, v in coeff_dict.items():
        if v != 0:
            assert k in ret_dict
            assert coeff_dict[k] == ret_dict[k]
        else:
            assert k not in ret_dict
    for k in ret_dict.keys():
        assert coeff_dict[k] == ret_dict[k]


def parse_ret(lines, ascode,unformat):
    from schubmult.perm_lib import uncode, permtrim

    ret_dict = {}
    for line in lines:
        try:
            k, v = line.split("  ")
        except ValueError:
            print(f"{line=}")
            continue
        try:
            v = unformat(v)
        except ValueError:
            continue
        ret_dict[(permtrim(literal_eval(k)) if not ascode else (uncode(literal_eval(k))))] = (
            unformat(v)
        )
    return ret_dict


base_dir = "schubmult_q"

json_files_data_args = load_json_test_names(base_dir)


@pytest.mark.parametrize("json_file", json_files_data_args)
def test_with_same_args_exec(capsys, json_file):
    from schubmult.utils import get_json, load_json_test_names
    from schubmult.utils.parsing import parse_coeff
    from schubmult.poly_lib import GeneratingSet
    from schubmult.perm_lib import permtrim, uncode

    args = get_json(f"{base_dir}/{json_file}")
    print(f"{json_file=} {args=} input_data")
    from schubmult.scripts.schubmult_q._script import main

    mult = args["mult"]  # noqa: F841
    mulstring = args["mulstring"]  # noqa: F841

    perms = args["perms"]

    ascode = args["ascode"]
    disp_mode = args["disp_mode"]
    pr = args["pr"]  # noqa: F841

    print(f"{args=}")
    print(f"{args['cmd_line']=}")

    from latex2sympy2_extended import latex2sympy
    from symengine import sympify

    unformat = {"basic": lambda v: parse_coeff(v), "latex": lambda v: parse_coeff(str(latex2sympy(v)))}

    ret_dict = main(args["cmd_line"])
    lines = capsys.readouterr()
    lines = str(lines.out).split("\n")

    if disp_mode != "raw":
        ret_dict = parse_ret(lines, ascode, unformat[disp_mode])
    elif ascode:
        ret_dict = {tuple(uncode(k)): v for k, v in ret_dict.items()}
    v_tuple = permtrim(perms[1]) if not ascode else uncode(perms[1])

    input_dict = {permtrim(perms[0]) if not ascode else uncode(perms[0]): 1}
    print(f"{v_tuple=} {input_dict=}")
    assert_dict_good(v_tuple, input_dict, ret_dict)


if __name__ == "__main__":
    test_with_same_args_exec()
