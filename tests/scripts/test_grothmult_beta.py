"""Beta specialization for all Grothendieck CLIs."""

import importlib
from ast import literal_eval

import pytest

from schubmult import Gx, Permutation
from schubmult.symbolic import expand, sympify
from schubmult.utils.parsing import parse_coeff


@pytest.mark.parametrize("script", ["grothmult_py", "grothmult_q", "grothmult_double", "grothmult_q_double"])
@pytest.mark.parametrize("value", ["0", "-1", "1/2", "t + 1"])
@pytest.mark.parametrize("display", ["raw", "basic"])
def test_beta_specialization(script, value, display, capsys):
    module = importlib.import_module(f"schubmult._scripts.{script}")
    argv = [script, "1", "3", "2", "-", "2", "1", "3"]
    original = module.main([*argv, "--display-mode", "raw"])
    capsys.readouterr()
    beta_value = sympify(value)
    expected = {p: sympify(c).subs(Gx._beta, beta_value) for p, c in original.items()}
    expected = {p: c for p, c in expected.items() if c != 0}
    assert Permutation([3, 2, 1]) in original

    result = module.main([*argv, f"--beta={value}", "--display-mode", display])
    output = capsys.readouterr().out
    if display == "basic":
        result = {}
        for line in output.splitlines():
            left, right = line.strip().split("  ", maxsplit=1)
            coefficient, perm = (left, right) if script == "grothmult_py" else (right, left)
            result[Permutation(literal_eval(perm))] = parse_coeff(coefficient)
    assert result.keys() == expected.keys()
    for perm, coefficient in result.items():
        assert expand(sympify(coefficient) - expected[perm]) == 0


@pytest.mark.parametrize("script", ["grothmult_py", "grothmult_q", "grothmult_double", "grothmult_q_double"])
def test_beta_passed_to_multiplication(script, monkeypatch):
    module = importlib.import_module(f"schubmult._scripts.{script}")
    seen = []

    def multiply(coeff_dict, perm, *args, **kwargs):
        seen.append(kwargs.get("beta", args[0] if args else None))
        return coeff_dict

    monkeypatch.setattr(module, script, multiply)
    module.main([script, "2", "1", "-", "2", "1", "--beta", "1/2", "--display-mode", "raw"])
    assert seen == [sympify("1/2")]


@pytest.mark.parametrize("script", ["grothmult_double", "grothmult_q_double"])
@pytest.mark.parametrize("value", ["0", "-1", "1/2", "t + 1"])
def test_beta_with_simplify(script, value):
    module = importlib.import_module(f"schubmult._scripts.{script}")
    argv = [script, "2", "1", "-", "2", "1", "--mixed-var", "--display-mode", "raw"]
    original = module.main(argv)
    result = module.main([*argv, "--simplify", f"--beta={value}"])
    from schubmult.symbolic import sympify_sympy

    for perm, coefficient in result.items():
        difference = sympify(coefficient) - sympify(original[perm]).subs(Gx._beta, sympify(value))
        assert sympify_sympy(difference).cancel() == 0


@pytest.mark.parametrize("value", ["0", "-1", "1/2", "t + 1"])
def test_beta_with_positive_display(value):
    from schubmult._scripts.grothmult_double import main
    from schubmult.symbolic import sympify_sympy

    argv = ["grothmult_double", "2", "1", "-", "2", "1", "--mixed-var", "--display-mode", "raw"]
    original = main(argv)
    result = main([*argv, "--display-positive", f"--beta={value}"])
    assert isinstance(result, dict)
    assert result.keys() == original.keys()
    for perm, coefficient in result.items():
        difference = sympify(coefficient) - sympify(original[perm]).subs(Gx._beta, sympify(value))
        assert sympify_sympy(difference).cancel() == 0


@pytest.mark.parametrize("script", ["grothmult_q", "grothmult_q_double"])
@pytest.mark.parametrize("value", ["0", "-1", "1/2", "t + 1"])
def test_beta_with_parabolic_projection(script, value):
    module = importlib.import_module(f"schubmult._scripts.{script}")
    argv = [script, "1", "3", "2", "4", "-", "2", "4", "1", "3", "--parabolic", "2", "2", "--display-mode", "raw"]
    original = module.main(argv)
    result = module.main([*argv, f"--beta={value}"])
    from schubmult.symbolic import sympify_sympy

    expected = {perm: sympify(sympify_sympy(coefficient).cancel()).subs(Gx._beta, sympify(value)) for perm, coefficient in original.items()}
    expected = {perm: coefficient for perm, coefficient in expected.items() if coefficient != 0}
    assert result.keys() == expected.keys()
    for perm, coefficient in result.items():
        difference = sympify(coefficient) - expected[perm]
        assert sympify_sympy(difference).cancel() == 0
