from schubmult import *

def compify(rc):
    if rc.perm.inv == 0:
        return {}
    last_desc = max(rc.perm.descents()) + 1
    new_rc, strip = rc.last_descent_strip()
    ret = {last_desc: strip}
    ret.update(compify(new_rc))
    return ret


if __name__ == "__main__":
    import sys

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    for perm in perms:
        decomp_set = set()
        for rc in RCGraph.all_rc_graphs(perm, len(perm.trimcode)):
            decomp = compify(rc)
            # trimmed_rc = rc
            # while trimmed_rc.perm.inv > 0:
            #     row = None
            #     for i in range(trimmed_rc.perm.inv):
            #         a, b = trimmed_rc.left_to_right_inversion(i)
            #         if b == a + 1:
            #             row, _ = trimmed_rc.left_to_right_inversion_coords(i)
            #             break
            #     to_add, trimmed_rc = trimmed_rc.vertical_cut(row)
            #     decomp.append(to_add)
            #decomp.append(trimmed_rc)
            #print(f"RC Graph:\n{rc}\ndecomposes to:")
            # pl = Plactic()
            # for desc, strip in sorted(decomp.items(),reverse=True):
            key = tuple(sorted(decomp.items()))
                #pl = pl.rs_insert(*tuple(reversed(strip)))
            assert key not in decomp_set
            decomp_set.add(key)