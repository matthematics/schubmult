from ast import literal_eval

import pytest

from schubmult.utils.test_utils import get_json, load_json_test_names


def assert_dict_good(v_tuple, input_dict, ret_dict, parabolic):
    from schubmult import schubmult_q, apply_peterson_woodward

    coeff_dict = schubmult_q(input_dict, v_tuple)

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

    for k, v in coeff_dict.items():
        if v != 0:
            assert k in ret_dict
            assert coeff_dict[k] == ret_dict[k]
        else:
            assert k not in ret_dict
    for k in ret_dict.keys():
        assert coeff_dict[k] == ret_dict[k]


def parse_ret(lines, ascode,unformat):
    from schubmult import uncode

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
        ret_dict[((literal_eval(k)) if not ascode else (uncode(literal_eval(k))))] = (
            unformat(v)
        )
    return ret_dict


base_dir = "schubmult_q"

json_files_data_args = load_json_test_names(base_dir)


@pytest.mark.parametrize("json_file", json_files_data_args)
def test_with_same_args_exec(capsys, json_file):
    from schubmult.utils.test_utils import get_json, load_json_test_names
    from schubmult.utils.parsing import parse_coeff
    from schubmult import uncode

    args = get_json(f"{base_dir}/{json_file}")
    print(f"{json_file=} {args=} input_data")
    from schubmult._scripts.schubmult_q import main

    mult = args["mult"]  # noqa: F841
    mulstring = args["mulstring"]  # noqa: F841

    perms = args["perms"]

    ascode = args["ascode"]
    disp_mode = args["disp_mode"]
    parabolic = [int(bob) for bob in args["parabolic"]]
    pr = args["pr"]  # noqa: F841

    print(f"{args=}")
    print(f"{args['cmd_line']=}")

    from latex2sympy2_extended import latex2sympy
    from schubmult.symbolic import sympify

    unformat = {"basic": lambda v: parse_coeff(v), "latex": lambda v: parse_coeff(str(latex2sympy(v)))}

    ret_dict = main(args["cmd_line"])
    lines = capsys.readouterr()
    lines = str(lines.out).split("\n")

    if disp_mode != "raw":
        ret_dict = parse_ret(lines, ascode, unformat[disp_mode])
    elif ascode:
        ret_dict = {tuple(uncode(k)): v for k, v in ret_dict.items()}
    v_tuple = (perms[1]) if not ascode else uncode(perms[1])

    input_dict = {(perms[0]) if not ascode else uncode(perms[0]): 1}
    print(f"{v_tuple=} {input_dict=}")
    assert_dict_good(v_tuple, input_dict, ret_dict, parabolic)


if __name__ == "__main__":
    test_with_same_args_exec()
