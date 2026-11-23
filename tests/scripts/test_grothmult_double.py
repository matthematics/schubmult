"""Script tests for ``grothmult_double``.

Each JSON case in ``tests/scripts/data/grothmult_double`` records a command line (generated
with the hidden ``-g`` flag). The test replays it through ``main``, parses the printed
``perm  coeff`` lines (or takes the raw dict in ``--display-mode raw``), and checks the
expansion by an independent polynomial identity: ``sum_w c_w G_w(x; y)`` must equal the
product of the input double Grothendieck polynomials, taken in ``y`` (``same``) or in
alternating alphabets ``y``, ``z``, ``y``, ... (``--mixed-var``), as polynomials in ``x``,
``y``, ``z`` and ``beta``. With ``--display-positive`` the printed coefficients are the
posified forms and are checked the same way.
"""

from ast import literal_eval

import pytest

from schubmult.utils.parsing import parse_coeff
from schubmult.utils.test_utils import get_json, load_json_test_names

base_dir = "grothmult_double"

json_files_data_args = load_json_test_names(base_dir)


def parse_ret(lines, ascode):
    from schubmult import Permutation, uncode

    ret_dict = {}
    for line in lines:
        try:
            k, v = line.strip().split("  ", maxsplit=1)
        except ValueError:
            continue
        try:
            key = literal_eval(k)
        except (ValueError, SyntaxError):
            continue
        perm = uncode(list(key)) if ascode else Permutation(list(key))
        ret_dict[perm] = parse_coeff(v)
    return ret_dict


def assert_expansion_good(perms, ret_dict, same):
    import sympy

    from schubmult import DGx, GeneratingSet
    from schubmult.symbolic import sympify, sympify_sympy

    var2 = GeneratingSet("y")
    var3 = var2 if same else GeneratingSet("z")

    lhs = sympify(1)
    for i, perm in enumerate(perms):
        lhs = lhs * DGx(perm, var2 if i % 2 == 0 else var3).expand()
    rhs = sum((sympify(v) * DGx(w).expand() for w, v in ret_dict.items()), sympify(0))
    # coefficients are rational in y, beta: clear denominators before comparing
    assert sympy.cancel(sympify_sympy(lhs - rhs)) == 0
    for v in ret_dict.values():
        assert sympy.cancel(sympify_sympy(v)) != 0


@pytest.mark.parametrize("json_file", json_files_data_args)
def test_with_same_args_exec(capsys, json_file):
    from schubmult import Permutation, uncode
    from schubmult._scripts.grothmult_double import main

    args = get_json(f"{base_dir}/{json_file}")
    ascode = args["ascode"]
    disp_mode = args["disp_mode"]
    same = args["same"]

    ret_dict = main(args["cmd_line"])
    out = capsys.readouterr().out

    if disp_mode == "raw":
        assert isinstance(ret_dict, dict)
        ret_dict = {Permutation(list(k)): v for k, v in ret_dict.items()}
    else:
        ret_dict = parse_ret(out.split("\n"), ascode)
    assert ret_dict

    perms = [uncode(p) if ascode else Permutation(p) for p in args["perms"]]
    assert_expansion_good(perms, ret_dict, same)


def test_display_positive_output_does_not_depend_on_jobs(capsys):
    from schubmult._scripts.grothmult_double import main

    cmd = ["grothmult_double", "2", "1", "3", "-", "1", "3", "2", "--mixed-var", "--display-positive"]
    outputs = []
    for jobs in ("1", "3"):
        main([*cmd, "--jobs", jobs])
        captured = capsys.readouterr()
        outputs.append((captured.out, captured.err))
    assert outputs[0][0]
    assert outputs[0] == outputs[1]


def test_jobs_defaults_to_one():
    from schubmult.utils.argparse import schub_argparse

    args, _ = schub_argparse(
        "grothmult_double",
        "test",
        ["2", "1", "3", "-", "1", "3", "2", "--display-positive"],
        yz=True,
        coprod=False,
    )
    assert args.jobs == 1


@pytest.mark.parametrize("jobs", ["0", "-1", "many"])
def test_jobs_must_be_a_positive_integer(capsys, jobs):
    from schubmult._scripts.grothmult_double import main

    with pytest.raises(SystemExit) as excinfo:
        main(["grothmult_double", "2", "1", "3", "-", "1", "3", "2", "--display-positive", "--jobs", jobs])
    assert excinfo.value.code == 2
    assert "--jobs" in capsys.readouterr().err
