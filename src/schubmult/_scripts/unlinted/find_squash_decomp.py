from schubmult.combinatorics.rc_graph import RCGraph
from schubmult.combinatorics.permutation import Permutation


def find_squash_decomp(rc):
    """Find decomposition (rc_sn, rc_grass) such that:
    - rc_sn.perm has len(rc_sn.perm) <= len(rc)
    - rc_grass is len(rc)-grassmannian (single descent at position len(rc)-1, zero-indexed)
    - rc_sn is of minimal size (fewest inversions)
    - rc_sn.resize(len(rc_grass)).squash_product(rc_grass) == rc

    Strategy: recursively right_zero_act() on rc until we reach
    rc_sn.disjoint_union(rc_grass), then extract the pieces.
    """
    rc = rc.normalize()
    n = len(rc)

    if len(rc.perm) - 1 <= n:
        # perm already fits within S_n, so rc_sn = rc and rc_grass = identity
        return rc, RCGraph([()]).resize(n)

    perm_size = len(rc.perm) - 1  # max moved element minus 1

    # BFS through right_zero_act to find a decomposable state
    stack = {rc}
    seen = set()
    best = None

    while stack:
        working_rc = stack.pop()
        if working_rc in seen:
            continue
        # Search bound from existing squash_decomp
        if len(working_rc) > n + perm_size:
            continue
        seen.add(working_rc)

        # Check if working_rc decomposes as a disjoint union rc_sn ⊔ rc_grass
        # coset_decomp with parabolic S_n: generators s_1,...,s_{n-1}
        min_cos, residue = working_rc.perm.coset_decomp(*list(range(1, n)))

        # For a valid decomposition:
        # - residue lives in S_n (len(residue) <= n)
        # - min_cos fixes {1,...,n} and has at most one descent
        if (len(residue) <= n
            and len(min_cos.descents()) <= 1
            and all(min_cos[i] == i + 1 for i in range(n))):

            # Extract rc_sn: the crosses in columns < n
            rc_sn = RCGraph([tuple([a for a in row if a < n]) for row in working_rc]).resize(n)

            # Extract rc_grass: everything not in rc_sn
            rc_grass = RCGraph([()]).resize(len(working_rc))
            for idx in range(working_rc.perm.inv):
                row, col = working_rc.left_to_right_inversion_coords(idx)
                if not rc_sn.has_element(row, col):
                    rc_grass = rc_grass.toggle_ref_at(row, col)

            rc_grass = rc_grass.vertical_cut(n)[0]

            if rc_grass.perm.inv == 0 or rc_grass.perm.descents() == {n - 1}:
                if len(rc_sn.perm) <= n:
                    # Check minimality: pick the one with smallest rc_sn
                    if best is None or rc_sn.perm.inv < best[0].perm.inv:
                        best = (rc_sn, rc_grass)

        stack.update(working_rc.right_zero_act())

    if best is None:
        raise ValueError(f"Failed to find squash decomposition for {rc}")

    rc_sn, rc_grass = best
    # Verify
    reconstructed = rc_sn.resize(len(rc_grass)).squash_product(rc_grass)
    assert reconstructed == rc, f"Reconstruction failed: {rc_sn}.resize({len(rc_grass)}).squash_product({rc_grass}) = {reconstructed} != {rc}"
    return rc_sn, rc_grass


if __name__ == "__main__":
    from schubmult import uncode

    for n in range(2, 6):
        for perm in Permutation.all_permutations(n):
            if perm.inv == 0:
                continue
            for rc in RCGraph.all_rc_graphs(perm, len(perm.trimcode)):
                try:
                    rc_sn, rc_grass = find_squash_decomp(rc)
                except ValueError as e:
                    print(f"rc = {rc}  SKIP: {e}")
                    print()
                    continue
                print(f"rc = {rc}")
                print(f"  rc_sn   = {rc_sn}  (perm = {rc_sn.perm})")
                print(f"  rc_grass = {rc_grass}  (perm = {rc_grass.perm})")
                print(f"  n-grass? descents = {rc_grass.perm.descents()}, n = {len(rc)}")
                print(f"  verified: rc_sn.resize({len(rc_grass)}).squash_product(rc_grass) == rc ? {rc_sn.resize(len(rc_grass)).squash_product(rc_grass) == rc}")
                print()
