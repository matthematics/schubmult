import pytest
from tests._tests import get_json
from ast import literal_eval


def assert_dict_good(v_tuple, input_dict, ret_dict, coprod, indices):
    from schubmult.schubmult_py import schubmult, schub_coprod

    if coprod:
        coeff_dict = schub_coprod(v_tuple, indices)
    else:
        coeff_dict = schubmult(input_dict, v_tuple)
    for k, v in coeff_dict.items():
        if v != 0:
            assert k in ret_dict
            assert coeff_dict[k] == ret_dict[k]
        else:
            assert k not in ret_dict
    for k in ret_dict.keys():
        assert coeff_dict[k] == ret_dict[k]


def parse_ret(lines, ascode, coprod):
    from schubmult.perm_lib import uncode, permtrim

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
            ret_dict[(tuple(literal_eval(k)) if not ascode else tuple(uncode(literal_eval(k))))] = (
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
                    k1 = tuple(permtrim(uncode(k1)))
                    k2 = tuple(permtrim(uncode(k2)))
                k = (k1, k2)
            except ValueError:
                print(f"{line=}")
                continue
            v = int(v)
            ret_dict[k] = v
    return ret_dict


base_dir = "script_tests/data/schubmult_py"

json_files_data_args = [
    f"{base_dir}/test_gen_py",
    f"{base_dir}/test_gen_py2",
    f"{base_dir}/test_gen_py_coprod",
]


@pytest.mark.parametrize("json_file", json_files_data_args)
def test_with_same_args_exec(capsys, json_file):
    from schubmult.perm_lib import uncode, permtrim

    args = get_json(json_file)
    print(f"{json_file=} {args=} input_data")
    from schubmult.schubmult_py._script import main

    mult = args["mult"]  # noqa: F841
    mulstring = args["mulstring"]  # noqa: F841

    perms = args["perms"]

    ascode = args["ascode"]
    coprod = args["coprod"]
    pr = args["pr"]  # noqa: F841

    print(f"{args=}")
    print(f"{args['cmd_line']=}")
    main(args["cmd_line"])
    lines = capsys.readouterr()
    lines = str(lines.out).split("\n")

    ret_dict = parse_ret(lines, ascode, coprod)
    v_tuple = (
        (tuple(perms[1]) if not ascode else tuple(uncode(perms[1])))
        if not coprod
        else (tuple(perms[0]) if not ascode else tuple(uncode(perms[0])))
    )
    input_dict = {tuple(permtrim(perms[0])) if not ascode else tuple(permtrim(uncode(perms[0]))): 1}
    indices = tuple(perms[1])
    print(f"{v_tuple=} {input_dict=} {indices=}")
    print("BOOB ASS")
    assert_dict_good(v_tuple, input_dict, ret_dict, coprod, indices)


if __name__ == "__main__":
    test_with_same_args_exec()
