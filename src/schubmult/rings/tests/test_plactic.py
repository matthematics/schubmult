from schubmult.rings.rc_graph import RCGraph
from schubmult.rings.plactic import Plactic, NilPlactic
from schubmult.perm_lib import Permutation


def test_ed_corresp():
    n = 6
    perms = Permutation.all_permutations(n)
    for perm in perms:
        rc_graphs = list(RCGraph.all_rc_graphs(perm))
        for rc in rc_graphs:
            P, Q = rc.edelman_greene()
            assert P.perm == ~perm, f"ED {P=} perm fail {perm=} {P.perm=}"
            Q2 = Plactic().rs_insert(*Q.row_word)
            assert Q == Q2, f"{Q.row_word=} != {Q2.row_word=}"
