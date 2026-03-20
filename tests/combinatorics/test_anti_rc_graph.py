from schubmult.combinatorics.anti_rc_graph import AntiRCGraph
from schubmult.combinatorics.permutation import Permutation
from schubmult.combinatorics.rc_graph import RCGraph


def test_anti_reading_order_accessors():
    anti = AntiRCGraph([(), (2,), (3, 1)])

    assert anti.anti_reduced_word == (3, 1, 2)
    assert anti.anti_compatible_sequence == (1, 1, 2)
    assert anti.as_reduced_anticompatible() == ((3, 1, 2), (1, 1, 2))


def test_from_reduced_anticompatible_roundtrip():
    word = (3, 1, 2)
    seq = (1, 1, 2)

    anti = AntiRCGraph.from_reduced_anticompatible(word, seq)

    assert anti.as_reduced_anticompatible() == (word, seq)
    assert anti.anti_is_valid


def test_anti_permutation_from_anti_word():
    word = (3, 1, 2)
    seq = (1, 1, 2)

    anti = AntiRCGraph.from_reduced_anticompatible(word, seq)

    assert anti.anti_permutation == Permutation.ref_product(*word)
    assert anti.perm == ~Permutation.ref_product(*word)


def test_perm_word_reads_anti_grid_order():
    anti = AntiRCGraph([(), (2,), (3, 1)])

    # right-to-left, top-to-bottom on anti grid
    assert anti.perm_word == (2, 3, 1)
    assert anti.reduced_word == (2, 3, 1)


def test_rcgraph_conversion_helpers():
    rc = RCGraph([(3, 1), (2,), ()])

    anti = AntiRCGraph.from_rc_graph(rc)
    rc2 = anti.to_rc_graph()

    assert tuple(anti) == ((), (2,), (3, 1))
    assert rc2 == rc


def test_anti_drawing_algorithmic_matrix_values():
    anti = AntiRCGraph([(), (2,), (3, 1)])

    matrix = [[anti[i, j] for j in range(anti.cols)] for i in range(anti.rows)]
    assert matrix == [[None, None, None], [2, None, None], [1, None, 3]]


def test_anti_display_name_in_repr():
    anti = AntiRCGraph([(2,), ()])

    assert str(anti).startswith("AntiRCGraph")


def test_anti_cols_not_truncated():
    anti = AntiRCGraph([(), (), (5, 1)])

    assert anti.cols == 5
    matrix = [[anti[i, j] for j in range(anti.cols)] for i in range(anti.rows)]
    assert matrix == [[None, None, None, None, None], [None, None, None, None, None], [1, None, None, None, 5]]


def test_anti_getitem_no_rcgraph_fallback():
    anti = AntiRCGraph([(), (2,), (5, 1)])

    assert anti[0] == ()
    assert anti[1] == (2,)
    assert anti[2] == (5, 1)
    assert anti[:] == ((), (2,), (5, 1))


def test_user_evacuation_example_display_shape():
    anti = AntiRCGraph([(4,), (3, 2), (1,)])

    matrix = [[anti[i, j] for j in range(anti.cols)] for i in range(anti.rows)]
    assert matrix == [[None, 4], [2, 3], [1, None]]

