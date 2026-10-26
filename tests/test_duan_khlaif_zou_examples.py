"""Inventory, convention, finite projection, and output tests for the paper examples."""

import itertools
import runpy
from collections import Counter
from pathlib import Path
from types import SimpleNamespace

import pytest

from schubmult import ChangeOfVariables, Permutation, QuantumDoubleGrothendieckRing
from schubmult.abc import y
from schubmult.symbolic import sympify_sympy
from schubmult.symbolic.poly.schub_poly import _vars

q = _vars.q_var


@pytest.fixture(scope="module")
def script():
    path = Path(__file__).resolve().parents[1] / "examples" / "duan_khlaif_zou_qk.py"
    return SimpleNamespace(**runpy.run_path(str(path)))


@pytest.fixture(scope="module")
def ring(script):
    return script.make_ring()


def example(script, identifier):
    return next(e for e in script.examples() if e.identifier == identifier)


def assert_coefficients(actual, expected):
    expected = {Permutation(w): c for w, c in expected.items()}
    assert actual.keys() == expected.keys()
    for w, c in expected.items():
        assert sympify_sympy(actual[w] - c).cancel() == 0


def test_inventory_and_parabolic_classes(script):
    inventory = script.examples()
    assert len(inventory) == 156
    assert len({e.identifier for e in inventory}) == 156
    assert Counter(e.flag.key for e in inventory) == {"fl3": 15, "fl4": 10, "fl134": 66, "fl124": 65}
    for e in inventory:
        for w in (e.left, e.right):
            assert sorted(w) == list(range(1, e.flag.n + 1))
            assert Permutation(w).min_coset_rep(*e.flag.parabolic) == Permutation(w)
    for key in ("fl3", "fl134", "fl124"):
        selected = [e for e in inventory if e.flag.key == key]
        flag = selected[0].flag
        classes = [w for w in itertools.permutations(range(1, flag.n + 1)) if w != tuple(range(1, flag.n + 1)) and Permutation(w).min_coset_rep(*flag.parabolic) == Permutation(w)]
        complete = set(itertools.combinations_with_replacement(sorted(classes), 2))
        actual = {tuple(sorted((e.left, e.right))) for e in selected}
        missing = {((3, 2, 1, 4), (3, 2, 1, 4))} if key == "fl124" else set()
        assert complete - actual == missing
        assert not actual - complete
    repeated = [e for e in inventory if e.occurrences > 1]
    assert len(repeated) == 1
    assert (repeated[0].left, repeated[0].right, repeated[0].page, repeated[0].occurrences) == ((4, 1, 2, 3), (3, 1, 4, 2), 42, 2)


def test_ring_uses_requested_change_of_variables(ring):
    assert isinstance(ring, QuantumDoubleGrothendieckRing)
    assert isinstance(ring.coeff_genset, ChangeOfVariables)
    assert ring.beta == -1
    assert ring.coeff_genset[1] == 1 - 1 / y[1]


def test_all_examples_compute_in_the_correct_finite_flag(script, ring):
    for e in script.examples():
        answer = script.compute(e, ring)
        assert answer, e.identifier
        allowed = {y[i] for i in range(1, e.flag.n + 1)} | {q[i] for i in range(1, len(e.flag.dimensions) + 1)}
        for w, c in answer.items():
            assert len(w) <= e.flag.n, e.identifier
            assert w == w.min_coset_rep(*e.flag.parabolic), e.identifier
            assert c.free_symbols <= allowed, e.identifier
            denominator = sympify_sympy(c).cancel().as_numer_denom()[1]
            assert not any(sympify_sympy(q[i]) in denominator.free_symbols for i in range(1, e.flag.n)), e.identifier


def test_fl3_quantum_square(script, ring):
    actual = script.compute(example(script, "fl3-01"), ring)
    assert_coefficients(actual, {
        (2, 1, 3): 1 - y[2] / y[1],
        (3, 1, 2): y[2] / y[1],
        (): q[1] * y[2] / y[1],
        (1, 3, 2): -q[1] * y[2] / y[1],
    })


def test_original_discrepancy_has_four_terms_and_no_quantum_parameters(script, ring):
    actual = script.compute(example(script, "fl4-04"), ring)
    assert_coefficients(actual, {
        (4, 3, 2, 1): y[2] * y[3] / y[1] ** 2,
        (3, 4, 2, 1): y[2] * (y[1] - y[3]) / y[1] ** 2,
        (4, 2, 3, 1): y[3] * (y[1] - y[2]) / y[1] ** 2,
        (3, 2, 4, 1): (y[1] - y[2]) * (y[1] - y[3]) / y[1] ** 2,
    })


def test_finite_projection_drops_stable_terms(script, ring):
    e = example(script, "fl4-02")
    stable = ring(e.left) * ring(e.right)
    assert any(len(w) > 4 for w in stable)
    # Test the computed finite image, not the paper's printed RHS, which
    # labels the q_3 term by 1423 rather than the computed 3214.
    assert_coefficients(script.compute(e, ring), {
        (3, 2, 4, 1): 1 - y[4] / y[1],
        (3, 2, 1, 4): q[3] * y[4] / y[1],
    })


@pytest.mark.parametrize("identifier,expected", [
    ("fl134-06", {
        (2, 3, 4, 1): 1 - y[4] / y[1],
        (2, 1, 3, 4): q[2] * y[4] / y[1],
    }),
    ("fl124-02", {
        (1, 4, 2, 3): 1 - y[4] / y[2],
        (2, 4, 1, 3): y[4] / y[2],
        (): q[2] * y[4] / y[2],
        (2, 1, 3, 4): -q[2] * y[4] / y[2],
    }),
])
def test_kato_projection_and_quantum_indices(script, ring, identifier, expected):
    assert_coefficients(script.compute(example(script, identifier), ring), expected)


@pytest.mark.parametrize("identifier", ["fl3-01", "fl4-04", "fl134-06", "fl124-02"])
def test_non_equivariant_specialization(script, ring, identifier):
    e = example(script, identifier)
    actual = script.compute(e, ring, non_equivariant=True)
    expected = script.compute(e, ring)
    expected = {tuple(w): sympify_sympy(c).subs({sympify_sympy(y[i]): 1 for i in range(1, e.flag.n + 1)}).cancel() for w, c in expected.items()}
    assert_coefficients(actual, {w: c for w, c in expected.items() if c != 0})
    assert not any(y[i] in c.free_symbols for c in actual.values() for i in range(1, e.flag.n + 1))


def test_list_does_not_compute(script, monkeypatch, capsys):
    def forbid_ring():
        pytest.fail("--list must not multiply")
    monkeypatch.setitem(script.main.__globals__, "make_ring", forbid_ring)
    assert script.main(["--list", "--section", "fl4"]) == 0
    out = capsys.readouterr()
    assert len(out.out.splitlines()) == 10
    assert "fl4-04" in out.out
    assert not out.err


@pytest.mark.parametrize("window,n,expected", [
    ((), 3, (1, 2, 3)),
    ((), 4, (1, 2, 3, 4)),
    ((2, 1), 4, (2, 1, 3, 4)),
    ((1, 3, 2), 4, (1, 3, 2, 4)),
])
def test_fixed_tails_and_identity_are_printed_correctly(script, window, n, expected):
    assert script.padded_permutation(Permutation(window), n) == expected
    assert script.padded_permutation(window, n) == expected
    assert script.basis(Permutation(window), n) == f"G~_({', '.join(map(str, expected))})(x;y)"


@pytest.mark.parametrize("tex", [False, True])
def test_cli_output_file_and_notation(script, tmp_path, capsys, tex):
    path = tmp_path / ("examples.tex" if tex else "examples.txt")
    args = ["--example", "fl4-04", "--output", str(path)]
    if tex:
        args.append("--tex")
    assert script.main(args) == 0
    out = capsys.readouterr()
    assert out.out == ""
    assert "[1/1] fl4-04" in out.err
    content = path.read_text()
    assert "fl4-04" in content
    assert "p.34" in content
    if tex:
        assert content.startswith(r"\documentclass")
        assert content.endswith("\\end{document}\n")
        assert r"\widetilde{\mathfrak{G}}" in content
        assert r"(x;y)" in content
        assert "z_i ->" not in content
    else:
        assert "G~_(4, 3, 2, 1)(x;y)" in content
        assert "y_2*y_3/y_1**2" in content


@pytest.mark.parametrize("args", [
    ["--example", "not-an-example"],
    ["--example", "fl4-04", "--section", "fl3"],
    ["--list", "--tex"],
    ["--list", "--non-equivariant"],
])
def test_invalid_cli_requests_are_explicit(script, args):
    with pytest.raises(SystemExit) as error:
        script.main(args)
    assert error.value.code == 2
