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


def test_parabolic_rejects_non_minimal_coset_rep(capsys):
    from schubmult._scripts.grothmult_q import main

    assert main(["grothmult_q", "2", "1", "3", "-", "2", "1", "3", "--parabolic", "2", "1"]) == 1
    assert "minimal coset representative" in capsys.readouterr().out
    assert main(["grothmult_q", "4", "1", "2", "3", "-", "2", "1", "3", "--parabolic", "1", "2"]) == 1


def test_parabolic_kato_sl3(capsys):
    """Kato (arXiv:1906.09343, Section 3), ``G = SL(3)``, ``P`` with ``W_P = <s_2>`` (blocks 1, 2), non-equivariant:
    ``O^{s_1} * O^{s_1} = O^{s_2 s_1}`` in ``QK(P^2)``; ``beta``-homogeneous form ``G_{213}^2 = G_{312}``.
    """
    from schubmult import Permutation
    from schubmult._scripts.grothmult_q import main

    ret = main(["grothmult_q", "2", "1", "3", "-", "2", "1", "3", "--parabolic", "1", "2", "--display-mode", "raw"])
    assert {Permutation(list(k)): v for k, v in ret.items()} == {Permutation([3, 1, 2]): 1}


def test_parabolic_kato_grassmannian_pieri():
    """``O_1 * O_{21} = O_{22} + q O_0 - q O_1`` in ``QK(Gr(2,4))`` (Buch--Mihalcea), via Kato's projection
    from ``QK(Fl_4)`` with blocks (2, 2); ``beta``-homogeneous with ``deg beta = -1``, ``deg q = 4``.
    """
    from schubmult import Permutation
    from schubmult._scripts.grothmult_q import main
    from schubmult.abc import beta
    from schubmult.symbolic import expand, sympify
    from schubmult.symbolic.poly.schub_poly import _vars

    q = _vars.q_var[1]
    ret = main(["grothmult_q", "1", "3", "2", "4", "-", "2", "4", "1", "3", "--parabolic", "2", "2", "--display-mode", "raw"])
    ret = {Permutation(list(k)): sympify(v) for k, v in ret.items()}
    expected = {Permutation([3, 4, 1, 2]): sympify(1), Permutation([1, 2, 3, 4]): q, Permutation([1, 3, 2, 4]): beta * q}
    assert set(ret) == set(expected)
    for w in expected:
        assert expand(ret[w] - expected[w]) == 0
