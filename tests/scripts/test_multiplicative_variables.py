"""Multiplicative character coordinates for the double Grothendieck CLIs."""

import importlib

import pytest

from schubmult import ChangeOfVariables, DoubleGrothendieckRing, Permutation, QuantumDoubleGrothendieckRing
from schubmult.abc import x, y, z
from schubmult.mult.groth_quantum_double import apply_kato
from schubmult.symbolic import sympify, sympify_sympy


def assert_coefficients(actual, expected):
    actual = {w: sympify_sympy(c).cancel() for w, c in actual.items()}
    expected = {w: sympify_sympy(c).cancel() for w, c in expected.items()}
    actual = {w: c for w, c in actual.items() if c != 0}
    expected = {w: c for w, c in expected.items() if c != 0}
    assert actual.keys() == expected.keys()
    assert all((actual[w] - c).cancel() == 0 for w, c in expected.items())


@pytest.mark.parametrize("script", ["grothmult_double", "grothmult_q_double"])
@pytest.mark.parametrize("mixed", [False, True])
@pytest.mark.parametrize("simplify", [False, True])
def test_matches_explicit_change_of_variables(script, mixed, simplify):
    module = importlib.import_module(f"schubmult._scripts.{script}")
    argv = [script, "2", "1", "-", "2", "1", "-", "1", "3", "2", "--display-mode", "raw"]
    if mixed:
        argv.append("--mixed-var")
    if simplify:
        argv.append("--simplify")
    actual = module.main([*argv, "--multiplicative-variables"])
    original = module.main([*argv, "--beta=-1"])
    # Independent simultaneous substitution of the conventional coordinates.
    substitution = {a[i]: 1 - 1 / a[i] for a in (y, z) for i in range(1, 10)}
    assert_coefficients(actual, {w: sympify(c).subs(substitution) for w, c in original.items()})

    cls = QuantumDoubleGrothendieckRing if script == "grothmult_q_double" else DoubleGrothendieckRing
    left = cls(x, ChangeOfVariables(z, y, lambda i, target: 1 - 1 / target[i]), beta=-1)
    right = cls(x, ChangeOfVariables(y, z, lambda i, target: 1 - 1 / target[i]), beta=-1) if mixed else left
    expected = left([2, 1]) * right([2, 1]) * right([1, 3, 2])
    assert_coefficients(actual, expected)


@pytest.mark.parametrize("script", ["grothmult_double", "grothmult_q_double"])
def test_multiplication_receives_changed_alphabets_and_beta(script, monkeypatch):
    module = importlib.import_module(f"schubmult._scripts.{script}")
    seen = []

    def multiply(coeffs, perm, left, right, beta):
        seen.append((left, right, beta))
        return coeffs

    monkeypatch.setattr(module, script, multiply)
    module.main([script, "2", "1", "-", "2", "1", "--multiplicative-variables", "--mixed-var", "--display-mode", "raw"])
    assert len(seen) == 1
    left, right, beta = seen[0]
    assert isinstance(left, ChangeOfVariables)
    assert isinstance(right, ChangeOfVariables)
    assert left[1] == 1 - 1 / y[1]
    assert right[1] == 1 - 1 / z[1]
    assert beta == -1


def test_quantum_parabolic_projection_and_code_input():
    from schubmult._scripts.grothmult_q_double import main

    # 1342 has code 011; this product uses the incidence-flag q_2.
    argv = ["grothmult_q_double", "0", "0", "1", "-", "0", "1", "1", "--code",
            "--multiplicative-variables", "--simplify", "--parabolic", "1", "2", "1", "--display-mode", "raw"]
    actual = main(argv)
    ring = QuantumDoubleGrothendieckRing(x, ChangeOfVariables(z, y, lambda i, t: 1 - 1 / t[i]), beta=-1)
    expected = apply_kato(ring([1, 2, 4, 3]) * ring([1, 3, 4, 2]), [2], n=4, beta=-1)
    assert_coefficients(actual, expected)
    Permutation.print_as_code = False


@pytest.mark.parametrize("script", ["grothmult_double", "grothmult_q_double"])
def test_explicit_beta_minus_one_is_allowed(script):
    module = importlib.import_module(f"schubmult._scripts.{script}")
    argv = [script, "2", "1", "-", "2", "1", "--multiplicative-variables", "--display-mode", "raw"]
    assert_coefficients(module.main([*argv, "--beta=-1"]), module.main(argv))


@pytest.mark.parametrize("script", ["grothmult_double", "grothmult_q_double"])
@pytest.mark.parametrize("options", [["--beta=0"], ["--beta=1/2"], ["--beta=t"], ["--display-positive", "--mixed-var"]])
def test_conflicting_options_fail_explicitly(script, options, capsys):
    module = importlib.import_module(f"schubmult._scripts.{script}")
    with pytest.raises(SystemExit) as exc:
        module.main([script, "2", "1", "-", "2", "1", "--multiplicative-variables", *options])
    assert exc.value.code == 2
    assert "--multiplicative-variables" in capsys.readouterr().err


@pytest.mark.parametrize("script", ["schubmult_double", "schubmult_q_double", "grothmult_py", "grothmult_q"])
def test_other_scripts_do_not_accept_the_flag(script, capsys):
    module = importlib.import_module(f"schubmult._scripts.{script}")
    with pytest.raises(SystemExit) as exc:
        module.main([script, "2", "1", "-", "2", "1", "--multiplicative-variables"])
    assert exc.value.code == 2
    assert "unrecognized arguments: --multiplicative-variables" in capsys.readouterr().err


@pytest.mark.parametrize("script", ["grothmult_double", "grothmult_q_double"])
@pytest.mark.parametrize("mode", ["basic", "sympy", "latex"])
def test_formatted_output_uses_multiplicative_characters(script, mode, capsys):
    module = importlib.import_module(f"schubmult._scripts.{script}")
    module.main([script, "2", "1", "-", "2", "1", "--multiplicative-variables", "--simplify", "--display-mode", mode])
    result = capsys.readouterr()
    assert result.out
    assert not result.err
    assert "beta" not in result.out and "β" not in result.out
    assert "_groth" not in result.out
