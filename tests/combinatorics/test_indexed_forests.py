from schubmult.combinatorics.indexed_forests import IndexedForest, Node, weak_composition_to_indfor


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


def test_trim_descent_directly_from_indexed_forest():
    forest = weak_composition_to_indfor((0, 2, 1, 0, 1, 0, 2))

    trimmed_at_5 = forest.trim_descent(5)
    # Inverse of blossoming at i=5: decrement c_5 and delete c_6.
    assert trimmed_at_5.code == (0, 2, 1, 0, 0, 2)
    assert trimmed_at_5.trim_descents == (3, 6)


def test_trim_descent_rejects_non_descent_index():
    forest = weak_composition_to_indfor((0, 2, 1, 0, 1, 0, 2))
    try:
        forest.trim_descent(4)
    except ValueError as exc:
        assert "not a trim descent" in str(exc)
    else:
        raise AssertionError("Expected ValueError for non-descent trim index")


def test_is_left_child_for_leaf_labels():
    #      2
    #     / \
    #    1   3
    n1 = Node(1)
    n2 = Node(2)
    n3 = Node(3)
    n2.left = n1
    n2.right = n3
    forest = IndexedForest((n2,))

    assert forest.is_left_child(1) is True
    assert forest.is_left_child(3) is False


def test_is_left_child_returns_false_for_non_leaf():
    #      3
    #     /
    #    2
    #   /
    #  1
    n1 = Node(1)
    n2 = Node(2)
    n3 = Node(3)
    n2.left = n1
    n3.left = n2
    forest = IndexedForest((n3,))

    # Node 2 is a left child of 3 but is not a leaf.
    assert forest.is_left_child(2) is False


def test_is_left_child_missing_label_raises():
    forest = weak_composition_to_indfor((1, 0, 0))
    try:
        forest.is_left_child(9)
    except ValueError as exc:
        assert "No node with index" in str(exc)
    else:
        raise AssertionError("Expected ValueError for missing node label")
