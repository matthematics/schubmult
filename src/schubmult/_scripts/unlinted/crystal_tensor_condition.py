from schubmult import RCGraph
from sympy import pretty_print


def _partitions(k, max_part):
    """Yield nonincreasing k-tuples of nonnegative ints summing to total.
    Parts are produced in nonincreasing order (standard partition convention)."""
    if max_part < 0:
        return
    if k < 0:
        return
    if k == 0:
        yield ()

    # choose next part v in range min(max_part, total) .. 0
    for v in range(max_part, -1, -1):
        for tail in _partitions(k - 1, v):
            yield (v,) + tail

def _fix(weight):
    """Given a weight (tuple of ints), return the list of descent positions."""
    descents = []
    for i in range(len(weight) - 1):
        if weight[i] == weight[i + 1]:
            descents.append(i + 1)
    if weight[-1] == 0:
        descents.append(len(weight))
    return descents

if __name__ == "__main__":
    import itertools
    import sys

    from schubmult import *
    n = int(sys.argv[1])
    m = int(sys.argv[2]) if len(sys.argv) > 2 else n

    perms = Permutation.all_permutations(n)
    weights = list(_partitions(m, m))
    bong = Permutation.w0(n - 1)
    for weight1 in weights:
        if sum(weight1) > 1 or sum(weight1) == 0:
            continue
        for perm in perms:
            for rc in RCGraph.all_rc_graphs(perm, length=n - 1):
                if not rc.is_highest_weight:
                    continue
                weight2 = rc.crystal_weight
                if sum(weight2) == 0:
                    continue
                rc_lw, _ = rc.to_lowest_weight()
                dem_perm = Permutation.sorting_perm(list(reversed(rc_lw.crystal_weight)))
                _coset_lambd = _fix(weight2)
                _coset_mu = _fix(weight1)
                vlam = dem_perm.min_coset_rep(*_coset_lambd)
                wmu = bong.max_coset_rep(*_coset_mu)
                if vlam.bruhat_leq(Permutation.longest_element(*(~wmu).descents(zero_indexed=False))):
                    #print(f"Weight {weight1} with Demazure weight {bong} and weight {weight2} with Demazure weight {dem_perm} does not pass the crystal tensor condition")
                    print("Good Pieri")
                    pretty_print(rc)
                else:
                    print("Bad")

