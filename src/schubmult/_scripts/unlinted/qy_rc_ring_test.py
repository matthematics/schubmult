from schubmult import *
from schubmult.rings.combinatorial import *
import itertools
import argparse

r = RCGraphRing()
qr = QYRCGraphRing()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("n", type=int)
    args = parser.parse_args()

    n = args.n

    # Generate all permutations of size n
    perms = Permutation.all_permutations(n)

    # Generate all RC graphs for each permutation
    # all_rc_graphs = [(perm, [RCGraph(perm, i) for i in range(len(perm.rc_graphs()))]) for perm in all_perms]

    # Loop over all pairs of permutations and pairs of RC graphs
    for perm1, perm2 in itertools.product(perms, repeat=2):
        for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1), RCGraph.all_rc_graphs(perm2)):
            elem1 = sum([coeff * qr._snap_qy(r(rc)) for rc, coeff in (r(rc1) * r(rc2)).items()])
            elem2 = qr._snap_qy(qr(rc1)) * qr._snap_qy(qr(rc2))
            assert elem1.almosteq(elem2), f"Failed for {perm1}, {perm2}, got elem1={elem1} and elem2={elem2}"