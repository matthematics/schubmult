if __name__ == "__main__":
    from schubmult import *
    from sympy import pretty_print
    import sys

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    rc_length = n - 1
    rc_ring = RCGraphRing()

    def raise_tensor_elem(rc_tensor_ring_elem, i):
        t_elem = rc_tensor_ring_elem.ring.zero
        for (rc1, rc2), coeff in rc_tensor_ring_elem.items():
            telm = CrystalGraphTensor(rc1, rc2).raising_operator(i)
            if telm is not None:
                t_elem += coeff * rc_tensor_ring_elem.ring((telm.factors[0], telm.factors[1]))
        return t_elem

    for perm in perms:
        for rc in RCGraph.all_rc_graphs(perm, rc_length):
            uplog = rc_ring(rc).coproduct()
            for i in range(1, rc.crystal_length()):
                rc_raised = rc.raising_operator(i)
                uplog_raised = uplog
                while rc_raised is not None:
                    uplog_raised = raise_tensor_elem(uplog_raised, i)
                    assert len(uplog) == len(uplog_raised)
                    print(f"I raised it {i}")
                    pretty_print(uplog_raised)
                    test = uplog - uplog_raised
                    assert any(v != 0 for v in test.values()), "Sanity check failed"
                    print("Good")
                    rc_raised = rc_raised.raising_operator(i)
            # hw = rc_.to_highest_weight()[0]
            # print("I got the highest weight")
            # pretty_print(hw)
            # assert len(hw) == len(uplog), "Highest weight length mismatch"