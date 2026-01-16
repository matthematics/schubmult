
def test_rc_bpd_ring_multiplication():
    from schubmult import RCGraph, BPD, Permutation, RCGraphRing, BPDRing
    from schubmult.symbolic import S
    import itertools
    n = 3
    perms = Permutation.all_permutations(n)
        
    ring1 = RCGraphRing()
    ring2 = BPDRing()
    for perm1, perm2 in itertools.product(perms, repeat=2):
        for len1 in range(len(perm1.trimcode), n):
            for len2 in range(len(perm2.trimcode), n):
                for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1, len1), RCGraph.all_rc_graphs(perm2, len2)):
                    rc_elem = ring1(rc1) * ring1(rc2)
                    bpd_elem = ring2(BPD.from_rc_graph(rc1)) * ring2(BPD.from_rc_graph(rc2))
                    assert all(v == S.Zero for v in (rc_elem - bpd_elem.to_rc_graph_ring_element()).values()), f"Error: RC graph ring element multiplication mismatch for permutations {perm1}, {perm2}:\nRC1:\n{rc1}\nRC2:\n{rc2}\nRC elem:\n{rc_elem}\nBPD elem:\n{bpd_elem.to_rc_graph_ring_element()}"
