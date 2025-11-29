"""
Check behaviour of RCGraph product vs. Edelman–Greene Q-tableau.

Usage:
    python -m schubmult.scripts.check_q_tableau N
This script is exploratory / best-effort.
"""
import itertools
import sys
from collections import Counter
from itertools import zip_longest

from schubmult import *


def main(argv):
    n = int(argv[1]) if len(argv) > 1 else 3

    perms = Permutation.all_permutations(n)
    for perm in perms:
        if perm.inv == 0:
            continue
        hw = set()
        for rc in RCGraph.all_rc_graphs(perm, length=n - 1):
            hw.add(rc.to_highest_weight()[0])
        bobset = {NilPlactic().ed_insert(*RootTableau.from_rc_graph(rc).reduced_word) for rc in hw}
        assert len(bobset) == len(hw)
    print("All checks passed")
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
