from schubmult.combinatorics.indexed_forests import weak_composition_to_indfor


def test_trim_descents_matches_qdes_rule_simple():
    # qdes(c) = {i : c_i > 0 and c_{i+1} = 0}
    forest = weak_composition_to_indfor((1, 1, 0, 0, 1, 0, 0))
    assert forest.trim_descents == (2, 5)
    assert tuple(node.index for node in forest.terminal_nodes) == forest.trim_descents
    assert tuple(node.index for node in forest.trim_descent_nodes) == forest.trim_descents


def test_trim_descents_matches_qdes_rule_with_gaps():
    forest = weak_composition_to_indfor((1, 0, 0, 1, 0, 1, 0, 0))
    assert forest.trim_descents == (1, 4, 6)


def test_trim_descents_matches_reported_nadeau_tewari_example():
    forest = weak_composition_to_indfor((0, 2, 1, 0, 1, 0, 2))
    assert forest.trim_descents == (3, 5, 7)
