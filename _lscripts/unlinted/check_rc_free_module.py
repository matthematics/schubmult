from schubmult import *

def partitions(n, k):
     if k == 0:
         return [[]]
     parts = []
     for p in partitions(n, k - 1):
         max_part = n if not p else p[-1]
         for i in range(min(max_part, n) + 1):
             parts.append(p + [i])
     return parts

def grass_perm_from_partition(partition, k):
    cd = [0] * (k - len(partition)) + list(reversed(partition))
    return uncode(cd)

def generate_grassmannian_permutations(n, k):
    grass_perms = set()
    parts = partitions(k, n)
    for part in parts:
        gp = grass_perm_from_partition(part, n)
        grass_perms.add(gp)
    return grass_perms

if __name__ == "__main__":
    import sys

    n = int(sys.argv[1])
    k = int(sys.argv[2])
    perms = Permutation.all_permutations(n)

    grass_perms = set()
    for pp in perms:
        if len(pp.descents()) == 1 and max(pp.descents(), default=-1) + 1 == n - 1:
            grass_perms.add(pp)
    seen_butts = set()
    butts_seen = {}
    for perm in perms:
        if len(perm) == n:
            continue
        for rc in RCGraph.all_rc_graphs(perm, n - 1):
            for grass_perm in grass_perms:
                # assert max(grass_perm.descents(), default=-1) + 1 == n or grass_perm.inv == 0, f"Not grassmannian: {grass_perm}"
                for grc in RCGraph.all_rc_graphs(grass_perm, n - 1):
                    new_rc = rc.squash_product(grc)
                    if len(new_rc.perm) > n + 1:
                        continue
                    try:
                        assert new_rc not in seen_butts, f"Duplicate RC graph found: {new_rc} {butts_seen[new_rc]} and {(rc, grc)}"
                    except AssertionError as e:
                        print(e)
                    butts_seen[new_rc] = butts_seen.get(new_rc, [])
                    butts_seen[new_rc].append((rc, grc))
                    seen_butts.add(new_rc)
    print("Bad ones:")
    for rc_graph in butts_seen:
        if len(butts_seen[rc_graph]) > 1:
            print(f"{rc_graph=}: {len(butts_seen[rc_graph])}")