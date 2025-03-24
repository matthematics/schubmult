import pytest
from schubmult._tests import get_json, load_json_test_names
from ast import literal_eval


def assert_dict_good(v_tuple, input_dict, ret_dict):
    from schubmult.schubmult_q import schubmult

    coeff_dict = schubmult(input_dict, v_tuple)
    for k, v in coeff_dict.items():
        if v != 0:
            assert k in ret_dict
            assert coeff_dict[k] == ret_dict[k]
        else:
            assert k not in ret_dict
    for k in ret_dict.keys():
        assert coeff_dict[k] == ret_dict[k]


def parse_ret(lines, ascode):
    from schubmult.perm_lib import uncode
    from symengine import sympify

    ret_dict = {}
    for line in lines:
        try:
            k, v = line.split("  ")
        except ValueError:
            print(f"{line=}")
            continue
        try:
            v = sympify(v)
        except ValueError:
            continue
        ret_dict[(tuple(literal_eval(k)) if not ascode else tuple(uncode(literal_eval(k))))] = (
            sympify(v)
        )
    return ret_dict


base_dir = "schubmult_q"

json_files_data_args = load_json_test_names(base_dir)


@pytest.mark.parametrize("json_file", json_files_data_args)
def test_with_same_args_exec(capsys, json_file):
    from schubmult.perm_lib import uncode, permtrim

    args = get_json(json_file)
    print(f"{json_file=} {args=} input_data")
    from schubmult.schubmult_q._script import main

    mult = args["mult"]  # noqa: F841
    mulstring = args["mulstring"]  # noqa: F841

    perms = args["perms"]

    ascode = args["ascode"]
    pr = args["pr"]  # noqa: F841

    print(f"{args=}")
    print(f"{args['cmd_line']=}")
    main(args["cmd_line"])
    lines = capsys.readouterr()
    lines = str(lines.out).split("\n")

    ret_dict = parse_ret(lines, ascode)
    v_tuple = (tuple(perms[1]) if not ascode else tuple(uncode(perms[1])))
        
    input_dict = {tuple(permtrim(perms[0])) if not ascode else tuple(permtrim(uncode(perms[0]))): 1}
    print(f"{v_tuple=} {input_dict=}")
    assert_dict_good(v_tuple, input_dict, ret_dict)


if __name__ == "__main__":
    test_with_same_args_exec()
