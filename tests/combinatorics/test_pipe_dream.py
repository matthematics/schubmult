from schubmult.combinatorics.permutation import Permutation
from schubmult.combinatorics.pipe_dream import PipeDream
from schubmult.combinatorics.rc_graph import RCGraph


def _first_rc(perm: Permutation, length: int):
    return next(iter(RCGraph.all_rc_graphs(perm, length=length)))


def test_pipe_dream_perm_word_simple_transposition_s1():
    rc = _first_rc(Permutation([2, 1]), length=2)
    pd = PipeDream.from_rc_graph(rc)
    assert pd.perm == rc.perm


def test_pipe_dream_perm_word_simple_transposition_s2():
    rc = _first_rc(Permutation([1, 3, 2]), length=3)
    pd = PipeDream.from_rc_graph(rc)
    assert pd.perm == rc.perm


def test_from_rc_graph_preserves_permutation_small_sizes():
    for n in range(1, 6):
        for perm in Permutation.all_permutations(n):
            for rc in RCGraph.all_rc_graphs(perm, length=n):
                assert PipeDream.from_rc_graph(rc).perm == rc.perm
