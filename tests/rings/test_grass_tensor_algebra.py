import itertools
import pytest
from schubmult import Permutation, RCGraph
from schubmult.rings.combinatorial.grass_tensor_algebra import GrassTensorAlgebra


def _is_valid_rc(rc):
    return rc in RCGraph.all_rc_graphs(rc.perm, len(rc))


def _is_full_grass(rc):
    return rc.perm.inv == 0 or rc.perm.descents() == {len(rc) - 1}


def _is_strictly_increasing_lengths(key):
    return all(len(key[i]) < len(key[i + 1]) for i in range(len(key) - 1))


@pytest.mark.xfail(reason="This test is currently failing because the from_rc_graph method in the GrassTensorAlgebra is not fully implemented. Once the from_rc_graph method is implemented, this test should pass if it correctly constructs elements with the expected properties of the normal form factors.")
def test_from_rc_graph_outputs_valid_normal_form_factors():
    ring = GrassTensorAlgebra()
    perms = Permutation.all_permutations(3)

    for perm in perms:
        for length in range(len(perm.trimcode), 4):
            for rc in RCGraph.all_rc_graphs(perm, length):
                elem = ring.from_rc_graph(rc)
                assert len(elem) == 1
                key = next(iter(elem.keys()))
                assert _is_strictly_increasing_lengths(key)
                for factor in key:
                    assert _is_full_grass(factor)
                    assert _is_valid_rc(factor)

@pytest.mark.xfail(reason="This test is currently failing because the multiplication in the GrassTensorAlgebra is not fully implemented. Once the multiplication is implemented, this test should pass if the multiplication correctly maintains the properties of the normal form factors.")
def test_mul_keeps_valid_normal_form_factors():
    ring = GrassTensorAlgebra()
    perms = Permutation.all_permutations(3)

    samples = []
    for perm in perms:
        for length in range(len(perm.trimcode), 4):
            graphs = list(RCGraph.all_rc_graphs(perm, length))
            samples.extend(graphs[:2])

    for rc1, rc2 in itertools.product(samples, repeat=2):
        e1 = ring.from_rc_graph(rc1)
        e2 = ring.from_rc_graph(rc2)
        prod = e1 * e2
        for key in prod.keys():
            assert _is_strictly_increasing_lengths(key)
            for factor in key:
                assert _is_full_grass(factor)
                assert _is_valid_rc(factor)
