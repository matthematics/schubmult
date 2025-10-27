from schubmult import *
from schubmult.rings.rc_graph import RCGraph
from schubmult.rings.rc_graph_ring import RCGraphRing
from schubmult.rings.crystal_graph import CrystalGraph, CrystalGraphTensor
import itertools


if __name__ == "__main__":
    from schubmult.schub_lib.schub_lib import elem_sym_perms
    from sympy import pretty_print
    import sys

    n = int(sys.argv[1])
    rc_ring = RCGraphRing()
    perms = Permutation.all_permutations(n)
    # for perm1, perm2 in itertools.product(perms, perms):
    #     graphs1 = RCGraph.all_rc_graphs(perm1, n-1)
    #     graphs2 = RCGraph.all_rc_graphs(perm2, n-1)
    #     dct = {}
    #     hw_check = set()
    #     for rc1, rc2 in itertools.product(graphs1, graphs2):
    #         cprd = CrystalGraphTensor(rc1, rc2).to_highest_weight()[0]
    #         if cprd in hw_check:
    #             continue
    #         hw_check.add(cprd)
    #         prod1 = rc_ring(rc1)*rc_ring(rc2)
    #         tring = rc_ring @ rc_ring
    #         raise_elem = tring.zero
    #         print(f"Coproduct of {perm1}, {perm2} highest weights")
    #         for tens in cprd.full_crystal:
    #             raise_elem = rc_ring.zero
    #             rc01, rc02 = tens.factors[0], tens.factors[1]
    #             prod2 = rc_ring(rc01)*rc_ring(rc02)
    #             for rc, coeff in prod2.items():
    #                 hw = CrystalGraphTensor(rc01, rc02).to_highest_weight()[0]
    #                 raise_elem += coeff * rc_ring(hw.factors[0]) * rc_ring(hw.factors[1])
    #             df = raise_elem-prod1
    #             print(f"{df=}")
    #             try:
    #                 assert all(v == 0 for v in df.values()), f"Mismatch in coproduct for {perm1}, {perm2} at highest weight {cprd} {raise_elem=} {prod1=}"

    class RCGraphCut(CrystalGraph):

        def __init__(self, under_rc, start, end):
            self.under_rc = under_rc
            self.start = start
            self.end = end

        def crystal_length(self):
            return self.end - self.start
        
        def epsilon(self, i):
            i += self.start
            if i < 1 or i >= self.end:
                return 0
            return self.under_rc.epsilon(i)
        
        def phi(self, i):
            i += self.start
            if i < 1 or i >= self.end:
                return 0
            return self.under_rc.phi(i)
        
        def raising_operator(self, index):
            index += self.start
            if index < 1 or index >= self.end:
                return None
            new_val = self.under_rc.raising_operator(index)
            if new_val is None:
                return None
            return RCGraphCut(new_val, self.start, self.end)

        def lowering_operator(self, index):
            index += self.start
            if index < 1 or index >= self.end:
                return None
            new_val = self.under_rc.lowering_operator(index)
            if new_val is None:
                return None
            return RCGraphCut(new_val, self.start, self.end)

        
            
    hw_elems = {}
    lr = {}
    tring = rc_ring @ rc_ring
    w0_rc = RCGraph.principal_rc(Permutation.w0(n), n-1)
    for perm in perms:
        
        #for rc in RCGraph.all_rc_graphs(perm, n-1):
        hw = set()
        results = set()
        for rc_p in RCGraph.all_rc_graphs(perm,n-1):
            prod = rc_ring(rc_p) * rc_ring(w0_rc)
            for rc, coeff in prod.items():
                rc_hw = CrystalGraphTensor(RCGraphCut(rc, 0, n-1), RCGraphCut(rc, n-1, len(rc))).reverse.to_highest_weight()[0].base_crystal
                full_rc_hw = RCGraph([*rc_hw.factors[0].under_rc.vertical_cut(n-1)[0], *rc_hw.factors[1].under_rc.rowrange(n-1)])
                if full_rc_hw in hw:
                    continue
            hw.add(full_rc_hw)
            results.add(full_rc_hw.vertical_cut(n-1)[0])
        for rc_result in results.keys():
            pretty_print(rc_result)
        input()
        # try:
        #     assert all(v == 0 for v in (cprd - our_cprd).values()), f"Mismatch in coproduct for perm {perm}"
        # except AssertionError as e:
        #     print(e)
        #     print(cprd - our_cprd)
        #     input()
        # for shift in range(len(prin_rc)):
        #     rc = RCGraph([*([()]* shift), *[tuple([a + len(prin_rc) for a in row]) for row in prin_rc], *([()]*(len(prin_rc) - shift))]).transpose().resize(2*len(prin_rc))
        #     i = len(prin_rc)
        #     hw = CrystalGraphTensor(RCGraphCut(rc, 0, i), RCGraphCut(rc, i, len(rc))).to_highest_weight()
        #     hw_rc = RCGraph([*hw[0].factors[0].under_rc.rowrange(0, i), *hw[0].factors[1].under_rc[i:]])
        #     rc1 = hw[0].factors[0].under_rc.vertical_cut(i)[0]
            
        #     rc2 = hw[0].factors[1].under_rc.rowrange(i)
        #     if (rc1, rc2) in the_set: 
        #         continue
        #     lr[perm] = lr.get(perm, tring.zero) + tring((rc1, rc2))
        #     the_set.add((rc1, rc2))
        #     hw_elems[(rc1, rc2)] = hw_elems.get((rc1, rc2), set())
        #     hw_elems[(rc1, rc2)].add(hw_rc)
            
        #     hw2 = CrystalGraphTensor(*rc.vertical_cut(i)).to_highest_weight()
        #     assert hw2[0].factors == (rc1, rc2)
    # for (rc1, rc2), elems in hw_elems.items():
    #     print(f"Highest weights for")
    #     pretty_print(CrystalGraphTensor(rc1, rc2))
    #     print(lr[(~rc1.perm, ~rc2.perm)])
    for perm in lr:
        elem = lr[perm]
        print(f"Perm {perm=}")
        pretty_print(lr[perm])
        print(ASx(~perm, n-1).coproduct())
            # for elem in elems:
            #     assert elem.is_valid
                
