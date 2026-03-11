from schubmult import *

if __name__ == "__main__":
    import itertools
    import sys
    from sympy import pretty_print
    
    # decompose n-row RC graph
    def decompose_rc(rc):
        """Decompose an n-row RC graph into a pair of (n-1)-row RC graph in S_n and a grass."""
        n = len(rc)
        stack = [rc]
        #def parabolic(k, perm):
        #rng = list(range(1, n))
        #    return perm.coset_decomp(*rng)
        seen = set()
        while stack:
            working_rc = stack.pop()
            # if len(working_rc) > 2 * n:
            #     print(stack)
            #     print(working_rc)
            #     print(seen)
            #     raise ValueError("should never get here")
            seen.add(working_rc.perm)
            #assert len(working_rc) == m
            min_cos, residue = working_rc.perm.coset_decomp(*list(range(1, len(working_rc))))
            if residue.inv == 0:
                return RCGraph([()]).resize(n), working_rc.vertical_cut(n)[0]
            if min_cos.inv == 0:
                return rc, RCGraph([()]).resize(n)
            #assert len(residue) <= n
            print(f"{residue=}, {min_cos=}")
            if len(residue) <= n and len(min_cos.descents()) <= 1 and min_cos * residue == residue * min_cos:
                working_rc2 = working_rc
                for a in reversed(residue.code_word):
                    working_rc2 = working_rc2.exchange_property(a)
                #ret_rc = working_rc.transpose().resize(len(residue.trimcode)).transpose().resize(n)
                ret_rc = RCGraph([tuple([a for a in row if a < n]) for row in working_rc]).resize(n)
                grass_rc = working_rc2
                if all(a > n for row in grass_rc for a in row):
                    #continue
                    grass_rc = grass_rc.vertical_cut(n)[0]
                    assert ret_rc.squash_product(grass_rc) == rc, "decomposition failed\n{}\n{}\n{}\n{}".format(rc, working_rc, ret_rc,grass_rc)
                    return ret_rc, grass_rc
            stack.extend([the_rc.normalize() for the_rc in working_rc.product(RCGraph([()])).keys() if the_rc.perm not in seen and len(the_rc.perm.coset_decomp(*list(range(1, len(the_rc) - 1)))[1]) <= n])
            seen.update([the_rc.perm for the_rc in stack])
            #pretty_print(stack)
        raise ValueError("should never get here")
    n = 6
    m = 4
    perms = [perm for perm in Permutation.all_permutations(n) if len(perm.trimcode) <= m]
    for perm in perms:
        for rc in RCGraph.all_rc_graphs(perm, m):
            print("rc:")
            pretty_print(rc)
            print("decomposes into:")
            pretty_print(decompose_rc(rc))
            print()
    