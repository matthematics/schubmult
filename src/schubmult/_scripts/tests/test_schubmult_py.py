from ast import literal_eval

import pytest

from schubmult.utils.test_utils import get_json, load_json_test_names


def assert_dict_good(v_tuple, input_dict, ret_dict, coprod, indices):
    import sys

    from schubmult import schub_coprod_py, schubmult_py
    if coprod:
        coeff_dict = schub_coprod_py(v_tuple, indices)
        print(f"{ret_dict=}",file=sys.stderr)
        print(f"{coeff_dict=}",file=sys.stderr)
    else:
        coeff_dict = schubmult_py(input_dict, v_tuple)
    for k, v in coeff_dict.items():
        if v != 0:
            assert coeff_dict[k] == ret_dict[k]
        else:
            assert k not in ret_dict
    for k in ret_dict.keys():
        assert coeff_dict[k] == ret_dict[k]


def parse_ret(lines, ascode, coprod):
    from schubmult import uncode

    ret_dict = {}
    if not coprod:
        for line in lines:
            try:
                v, k = line.split("  ")
            except ValueError:
                print(f"{line=}")
                continue
            try:
                v = int(v)
            except ValueError:
                continue
            ret_dict[((literal_eval(k)) if not ascode else uncode(literal_eval(k)))] = (
                int(v)
            )
    else:
        for line in lines:
            line = str(line)
            first_split = ") (" if not ascode else "] ["
            second_split = " (" if not ascode else " ["
            jn = "),(" if not ascode else "],["
            try:
                vf, s = line.split(first_split)
                v, f = vf.split(second_split)
                k1, k2 = literal_eval(f"({second_split}{f + jn + s})")
                if ascode:
                    k1 = ((uncode(k1)))
                    k2 = ((uncode(k2)))
                k = (k1, k2)
            except ValueError:
                print(f"{line=}")
                continue
            v = int(v)
            ret_dict[k] = v
    return ret_dict


base_dir = "schubmult_py"

json_files_data_args = load_json_test_names(base_dir)


@pytest.mark.parametrize("json_file", json_files_data_args)
def test_with_same_args_exec(capsys, json_file):
    from schubmult import uncode
    from schubmult import Permutation

    args = get_json(f"{base_dir}/{json_file}")
    print(f"{json_file=} {args=} input_data")
    from schubmult._scripts.schubmult_py import main

    mult = args["mult"]  # noqa: F841
    mulstring = args["mulstring"]  # noqa: F841

    perms = args["perms"]
    disp_mode = args["disp_mode"]
    ascode = args["ascode"]
    coprod = args["coprod"]
    pr = args["pr"]  # noqa: F841

    print(f"{args=}")
    print(f"{args['cmd_line']=}")
    ret_dict = main(args["cmd_line"])
    lines = capsys.readouterr()
    lines = str(lines.out).split("\n")

    if disp_mode != "raw":
        ret_dict = parse_ret(lines, ascode, coprod)
    elif coprod:
        if ascode:
            ret_dict = {((uncode(list(k[0]))),(uncode(list(k[1])))): v for k, v in ret_dict.items()}
    elif ascode:
        ret_dict = {(uncode(list(k))): v for k, v in ret_dict.items()}
    v_tuple = (
        (Permutation(perms[1]) if not ascode else (uncode(perms[1])))
        if not coprod
        else (Permutation(perms[0]) if not ascode else (uncode(perms[0])))
    )
    input_dict = {((perms[0])) if not ascode else ((uncode(perms[0]))): 1}
    indices = tuple(perms[1])
    print(f"{v_tuple=} {input_dict=} {indices=}")
    assert_dict_good(v_tuple, input_dict, ret_dict, coprod, indices)


if __name__ == "__main__":
    test_with_same_args_exec()
