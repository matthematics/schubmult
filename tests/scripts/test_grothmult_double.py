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
