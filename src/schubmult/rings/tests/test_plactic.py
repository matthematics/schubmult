from schubmult.rings.rc_graph import RCGraph
from schubmult.schub_lib.plactic import Plactic, NilPlactic
from schubmult.schub_lib.perm_lib import Permutation


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

import pytest

from schubmult.schub_lib.perm_lib import Permutation
from schubmult.schub_lib.plactic import Plactic, NilPlactic
from schubmult.rings.rc_graph import RCGraph


def check_one_rc(rc):
    # pick reasonable i range from entries
    pl = rc.weight_tableau
    entries = [e for row in pl._word for e in row]
    if not entries:
        return True
    max_e = max(entries)
    for i in range(1, max_e):
        # Plactic ops
        p_e = pl.raising_operator(i)
        p_f = pl.lowering_operator(i)

        # RCGraph ops via highest-weight RC
        rc_e = rc.raising_operator(i)
        rc_f = rc.lowering_operator(i)

        # compare None cases
        if (rc_e is None) != (p_e is None):
            return False, f"Mismatch e_i None state for i={i}: pl={p_e} rc={rc_e}"
        if (rc_f is None) != (p_f is None):
            return False, f"Mismatch f_i None state for i={i}: pl={p_f} rc={rc_f}"

        if rc_e is not None:
            if rc_e.p_tableau != p_e:
                return False, f"e_i mismatch for i={i}: pl={p_e} rc={rc_e.p_tableau}"
            # involutive check
            back = p_e.lowering_operator(i)
            if back != pl:
                return False, f"e then f failed for i={i}: got {back}, want {pl}"

        if rc_f is not None:
            if rc_f.p_tableau != p_f:
                return False, f"f_i mismatch for i={i}: pl={p_f} rc={rc_f.p_tableau}"
            back = p_f.raising_operator(i)
            if back != pl:
                return False, f"f then e failed for i={i}: got {back}, want {pl}"

        # shape preserved
        if p_e is not None and p_e.shape != pl.shape:
            return False, f"shape changed by e_i for i={i}: {pl.shape} -> {p_e.shape}"
        if p_f is not None and p_f.shape != pl.shape:
            return False, f"shape changed by f_i for i={i}: {pl.shape} -> {p_f.shape}"

    return True, ""


@pytest.mark.parametrize("n", [3, 4])
def test_plactic_vs_rcgraph_small(n):
    perms = Permutation.all_permutations(n)
    # iterate several RCGraphs (principal) to get a variety of tableaux
    for p in perms:
        for rc in RCGraph.all_rc_graphs(p, length=n - 1):
            ok = check_one_rc(rc)
            assert ok, f"Failure for rc {rc}"