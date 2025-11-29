if __name__ == "__main__":
    from functools import cache
    
    @cache
    def eg_weak_order(rc1, rc2):
        if rc1.perm.inv > rc2.perm.inv:
            return False
        if rc1.perm.inv == rc2.perm.inv:
            return RootTableau.from_rc_graph(rc1).edelman_greene_invariant == RootTableau.from_rc_graph(rc2).edelman_greene_invariant
        for d in rc2.perm.descents():
            real_descent = d + 1
            if eg_weak_order(rc1, rc2.exchange_property(real_descent)):
                return True
            for i in range(1, rc2.crystal_length()):
                rc2_try = rc2.lowering_operator(i)
                while rc2_try is not None:
                    if eg_weak_order(rc1, rc2_try):
                        return True
                    rc2_try = rc2_try.lowering_operator(i)
        return False


    from schubmult import *
    from sympy import pretty_print
    import sys

    n = int(sys.argv[1])

    # find algebraic presentation of RC graphs

    perms = Permutation.all_permutations(n)
    w0 = Permutation.w0(n)
    rc_length = n - 1
    w0_rc = RCGraph.principal_rc(w0, rc_length)
    for perm in perms:
        for rc in RCGraph.all_rc_graphs(perm, rc_length):
            try:
                assert eg_weak_order(rc, w0_rc)
            except AssertionError:
                print("bad")
                pretty_print(rc)
            print("good")

    