"""Script tests for ``grothmult_q_double``.

Each JSON case in ``tests/scripts/data/grothmult_q_double`` records a command line (generated
with the hidden ``-g`` flag).  The test replays it through ``main``, parses the printed
``perm  coeff`` lines (or takes the raw dict in ``--display-mode raw``), and checks the
expansion by an independent polynomial identity: ``sum_w c_w G^q_w(x; y)`` must equal the
product of the input quantum double Grothendieck polynomials, taken in ``y`` (``same``) or in
alternating alphabets ``y``, ``z``, ``y``, ... (``--mixed-var``), as polynomials in ``x``, ``y``,
``z``, ``q`` and ``beta``, with ``G^q_w = lm_quantize(G_w)`` (``qgroth_poly``).
"""

from ast import literal_eval

import pytest

from schubmult.utils.parsing import parse_coeff
from schubmult.utils.test_utils import get_json, load_json_test_names

base_dir = "grothmult_q_double"

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

    from schubmult import GeneratingSet
    from schubmult.abc import x
    from schubmult.mult.groth_quantum_double import qgroth_poly
    from schubmult.symbolic import sympify, sympify_sympy

    var2 = GeneratingSet("y")
    var3 = var2 if same else GeneratingSet("z")

    lhs = sympify(1)
    for i, perm in enumerate(perms):
        lhs = lhs * qgroth_poly(perm, x, var2 if i % 2 == 0 else var3)
    rhs = sum((sympify(v) * qgroth_poly(w, x, var2) for w, v in ret_dict.items()), sympify(0))
    # coefficients are rational in y, beta: clear denominators before comparing
    assert sympy.cancel(sympify_sympy(lhs - rhs)) == 0
    for v in ret_dict.values():
        assert sympy.cancel(sympify_sympy(v)) != 0


@pytest.mark.parametrize("json_file", json_files_data_args)
def test_with_same_args_exec(capsys, json_file):
    from schubmult import Permutation, uncode
    from schubmult._scripts.grothmult_q_double import main

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


@pytest.mark.parametrize("flag", [["--display-positive", "--mixed-var"], ["--parabolic", "2"], ["--nil-hecke", "2"]])
def test_unsupported_options_refused(capsys, flag):
    from schubmult._scripts.grothmult_q_double import main

    assert main(["grothmult_q_double", "2", "1", "-", "2", "1", *flag]) == 1
    assert "not supported" in capsys.readouterr().out
