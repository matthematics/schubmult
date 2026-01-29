from schubmult import *

def strip_replacement(rc: RCGraph) -> RCGraph:
    maxd = len(rc.perm.trimcode)
    if maxd == 0:
        return rc
    rows = []
    interim = rc
    while len(interim.perm.trimcode) >= maxd:
        interim, row = interim.exchange_property(max(interim.perm.descents(zero_indexed=False)), return_row=True)
        rows.append(row)
    interim = interim.kogan_kumar_insert(maxd, rows)
    return interim

if __name__ == "__main__":
    import sys

    n = int(sys.argv[1])
    
    for w in Permutation.all_permutations(n):
        for rc in RCGraph.all_rc_graphs(w, n):
            assert rc == strip_replacement(rc)