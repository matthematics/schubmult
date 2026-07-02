from schubmult.combinatorics.indexed_forests import (
    IndexedForest,
    Node,
    indexed_forest_from_trimming_word,
    letterpair,
    omega_setvalued_insertion,
    omega_setvalued_insertion_from_wcgraph,
    omega_is_set_valued_compatible,
    omega_set_valued_compatibility,
    weak_composition_to_indfor,
)


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


def test_indexed_forest_from_trimming_word_round_trip_example():
    trimming_word = (1, 1, 2, 4, 7)
    forest = indexed_forest_from_trimming_word(trimming_word)
    expected = weak_composition_to_indfor((2, 1, 0, 1, 0, 0, 1, 0))
    assert forest == expected


def test_indexed_forest_from_trimming_word_empty():
    forest = indexed_forest_from_trimming_word(())
    assert len(forest.roots) == 0
    assert forest.code == ()


def test_indexed_forest_from_trimming_word_rejects_nonpositive_entries():
    try:
        indexed_forest_from_trimming_word((1, 0, 2))
    except ValueError as exc:
        assert "positive integers" in str(exc)
    else:
        raise AssertionError("Expected ValueError for non-positive trimming-word entry")


def test_omega_set_valued_compatibility_true_case():
    # Longer compatible sequence reducing to one reduced letter with set {1,2}.
    pairs = (letterpair(1, 1), letterpair(1, 2))
    data = omega_set_valued_compatibility(pairs)
    assert data["ok"] is True
    assert data["omega_reduced_word"] == (1, 1)
    assert data["wc_reduced_word"] == (1,)
    assert data["set_sequence"] == ((1, 2),)
    assert omega_is_set_valued_compatible(pairs) is True


def test_omega_set_valued_compatibility_false_case():
    # Not weakly increasing in the compatible sequence, so WCGraph rejects it.
    pairs = (letterpair(1, 2), letterpair(1, 1))
    data = omega_set_valued_compatibility(pairs)
    assert data["ok"] is False
    assert omega_is_set_valued_compatible(pairs) is False


def test_omega_setvalued_insertion_merges_on_unreduced_step():
    # Word (1,1) is unreduced at second step; both compatible labels merge.
    data = omega_setvalued_insertion((1, 1), (1, 2))

    assert data["reduced_word"] == (1,)
    assert data["set_sequence"] == ((1, 2),)
    assert len(data["node_sets"]) == 1
    assert next(iter(data["node_sets"].values())) == (1, 2)


def test_omega_setvalued_insertion_from_wcgraph_matches_wc_reduction():
    from schubmult.combinatorics.wc_graph import WCGraph

    wc = WCGraph.from_word_compatible((1, 1), (1, 2), length=2)
    data = omega_setvalued_insertion_from_wcgraph(wc)

    red_word, set_seq = wc.to_reduced_compatible_set_sequence()
    assert data["reduced_word"] == red_word
    assert data["set_sequence"] == set_seq
