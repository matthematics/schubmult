"""Script tests for ``grothmult_q``.

Each JSON case in ``tests/scripts/data/grothmult_q`` records a command line (generated with
the hidden ``-g`` flag).  The test replays it through ``main``, parses the printed
``perm  coeff`` lines (or takes the raw dict in ``--display-mode raw``), and checks the
expansion by an independent polynomial identity: ``sum_w c_w G^q_w(x)`` must equal the
product of the input quantum Grothendieck polynomials as polynomials in ``x``, ``beta`` and
``q``, where ``G^q_w = lm_quantize(G_w)`` is the Lenart--Maeno quantization.
"""

from ast import literal_eval

import pytest

from schubmult.utils.parsing import parse_coeff
from schubmult.utils.test_utils import get_json, load_json_test_names

base_dir = "grothmult_q"

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


def assert_expansion_good(perms, ret_dict):
    from schubmult import Gx
    from schubmult.abc import x
    from schubmult.mult.groth_quantum_double import lm_quantize
    from schubmult.symbolic import expand, sympify

    beta = Gx._beta

    def gq(w):
        return lm_quantize(Gx(w).expand(), max(len(w), 2), x, beta)

    lhs = sympify(1)
    for perm in perms:
        lhs = lhs * gq(perm)
    rhs = sum((sympify(v) * gq(w) for w, v in ret_dict.items()), sympify(0))
    assert expand(lhs - rhs) == 0
    for v in ret_dict.values():
        assert expand(v) != 0


@pytest.mark.parametrize("json_file", json_files_data_args)
def test_with_same_args_exec(capsys, json_file):
    from schubmult import Permutation, uncode
    from schubmult._scripts.grothmult_q import main

    args = get_json(f"{base_dir}/{json_file}")
    ascode = args["ascode"]
    disp_mode = args["disp_mode"]

    ret_dict = main(args["cmd_line"])
    out = capsys.readouterr().out

    if disp_mode == "raw":
        assert isinstance(ret_dict, dict)
        ret_dict = {Permutation(list(k)): v for k, v in ret_dict.items()}
    else:
        ret_dict = parse_ret(out.split("\n"), ascode)
    assert ret_dict

    perms = [uncode(p) if ascode else Permutation(p) for p in args["perms"]]
    assert_expansion_good(perms, ret_dict)


def test_parabolic_not_supported(capsys):
    from schubmult._scripts.grothmult_q import main

    assert main(["grothmult_q", "2", "1", "-", "2", "1", "--parabolic", "2"]) == 1
    assert "not supported" in capsys.readouterr().out
