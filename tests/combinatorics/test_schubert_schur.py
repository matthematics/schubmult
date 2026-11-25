import pytest

from schubmult.combinatorics.experimental.schubert_schur import coefficient, diagram_count, product
from schubmult.combinatorics.permutation import Permutation, uncode
from schubmult.mult.single import schubmult_py


@pytest.mark.parametrize(
    ("lam", "w"),
    [
        ((3, 2, 1), (3, 4, 6, 7, 5, 1, 2)),  # Neither the shape nor its complement is a hook.
        ((3, 3, 2), (4, 5, 6, 7, 2, 3, 1)),  # A hook complement on a disconnected forest.
    ],
)
def test_mixed_descent_multiplicity(lam, w):
    u = (2, 1, 5, 3, 7, 6, 4)
    expected = schubmult_py({Permutation(u): 1}, uncode([0, *reversed(lam)]))[Permutation(w)]
    assert coefficient(u, lam, 4, w) == expected == 2


def test_entire_stable_product():
    u = (2, 1, 3, 4, 6, 5)
    expected = {w: c for w, c in schubmult_py({Permutation(u): 1}, uncode([0, 2, 3, 3])).items() if c}
    assert len(expected) == 7 and any(len(w) == 8 for w in expected)
    assert product(u, (3, 3, 2), 4) == expected
    assert product(u, (), 4) == {Permutation(u): 1}


def test_hook_cycle_and_empty_interval():
    square = [(1, 1), (1, 2), (2, 1), (2, 2)]
    assert diagram_count(square, (3, 1)) == 0
    assert diagram_count(square, (2, 2)) == 1
    assert coefficient((2, 1), (1,), 2, (1, 3, 4, 2)) == 0
