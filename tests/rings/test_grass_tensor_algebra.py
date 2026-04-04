import itertools
import pytest
from schubmult import Permutation, RCGraph, Sx
from schubmult.rings.combinatorial.grass_tensor_algebra import GrassTensorAlgebra


def test_match_rc_graph():
    ring = GrassTensorAlgebra()
    perms = Permutation.all_permutations(5)

    for perm in perms:
        cem = RCGraph.full_CEM(perm, len(perm.trimcode))
        for rc, cem_dict in cem.items():
            elem = ring.from_dict(cem_dict)
            rc_elem = elem.to_rc_graph_ring_element()
            if rc.perm == perm:
                assert rc_elem.almosteq(rc_elem.ring(rc)), f"Failed on {rc}\nGot {rc_elem}\nExpected {rc_elem.ring(rc)}"
            else:
                assert rc_elem.almosteq(rc_elem.ring.zero), f"Failed on {rc}\nGot {rc_elem}\nExpected 0"

def test_mul_matches_schubert():
    ring = GrassTensorAlgebra()
    n = 4
    perms = Permutation.all_permutations(n)

    grass_tensor_elems = {}
    for perm in perms:
        grass_tensor_elems[perm] = ring.zero
        cem = RCGraph.full_CEM(perm, n - 1)
        for rc, cem_dict in cem.items():
            grass_tensor_elems[perm] += ring.from_dict(cem_dict)
    for perm1, perm2 in itertools.product(perms, repeat=2):
        prd = Sx(perm1) * Sx(perm2)
        prd_elem = grass_tensor_elems[perm1] * grass_tensor_elems[perm2]
        prd_rc = prd_elem.to_rc_graph_ring_element()
        for rc, coeff in prd_rc.items():
            assert prd.get(rc.perm, 0) == coeff, f"Failed on {perm1} * {perm2}, got {rc.perm}: {coeff} which is not in {prd}\n{prd.get(rc.perm,0)=}\n{rc.perm=}"