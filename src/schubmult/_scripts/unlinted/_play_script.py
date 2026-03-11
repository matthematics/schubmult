from schubmult import *

if __name__ == "__main__":
    import itertools
    import sys
    from sympy import pretty_print
    
    # decompose n-row RC graph
    def bpd_squash(rc1, rc2):
        # bpd1, bpd2 = BPD.from_rc_graph(rc1.transpose()), BPD.from_rc_graph(rc2.transpose())
        # bpd1 = bpd1.resize(len(rc1))
        # bpd2 = bpd2.resize(len(rc2))
        # rc11 = bpd1.transpose().to_rc_graph().resize(len(rc1))
        # rc22 = bpd2.transpose().to_rc_graph().resize(len(rc2))

        squash = rc1.squash_product(rc2)
        # if rc11 != rc1 or rc22 != rc2:
        #     print("hoopy")
        #return rc11.squash_product(rc22)
        return BPD.from_rc_graph(squash).transpose().to_rc_graph().transpose().resize(len(rc1))

    def bpd_disjoint_union(bpd1, bpd2):
        pass

    
    def recurse_bpd(w, length, n):
        import numpy as np
        if length == 0:
            if w.inv == 0:
                return {BPD(np.array([], dtype=TileType))}
            return set()
        if length == 1:
            if len(w.trimcode) <= 1:
                arr = np.array([[TileType.BLANK if i < w.inv else (TileType.ELBOW_SE if i == w.inv else TileType.HORIZ) for i in range(len(w))]])
                return {BPD(arr)}
            return set()
        last = pull_out_var(length, w)
        ret = set()
        for num_set, perm2 in last:
            old_set = recurse_bpd(perm2, length - 1, n)
            for old_bpd in old_set:
                # recurse left
                new_grid = old_bpd._grid.copy()
                the_len = max(length, len(w))
                if new_grid.shape[1] < the_len:
                    extension = the_len - new_grid.shape[1]
                else:
                    extension = 0
                
                new_grid = np.pad(new_grid, ((0, 0), (0, extension)), constant_values=TileType.TBD)
                new_row = np.full(the_len, dtype=TileType, fill_value=TileType.TBD)
                new_grid = np.vstack([new_grid, new_row])
                perm = w * Permutation.w0(n - length).shiftup(length)
                new_grid[-1, :] = BPD.row_from_k_chain(perm2, perm, length, len(w))
                # for i in w[length:]:
                #     if w[i] == perm2[i - 1]:
                #         new_grid[-1, w[i] - 1] = TileType.CROSS
                # #new_grid[-1, :] = new_row
                # for j in range(len(new_grid)):
                #     if new_grid[j, -1] == TileType.VERT:
                #         if j == 0:
                #             if TileType(new_grid[j, -2]).feeds_right:
                #                 new_grid[j, -1] = TileType.HORIZ
                #             else:
                #                 new_grid[j, -1] = TileType.ELBOW_SE
                #         else:
                #             if TileType(new_grid[j, -2]).feeds_right:
                #                 if TileType(new_grid[j - 1, -1]).entrance_from_bottom:
                #                     new_grid[j, -1] = TileType.CROSS
                #                 else:
                #                     new_grid[j, -1] = TileType.HORIZ
                #             else:
                #                 new_grid[j, -1] = TileType.ELBOW_SE
                ret.add(BPD(new_grid))
        return ret

    def to_bruhat_path(self):
        n = len(self.perm)
        bigself = self.resize(n)
        return tuple([bigself.resize(i).perm * Permutation.w0(n - i).shiftup(i) for i in range(n - 1, -1, -1)])

    def bruhat_bpd(w):
        bpath = []
        w2 = w
        ret = set()
        stack = [(w2, bpath, len(w))]
        while stack:
            w3, bpath, i = stack.pop()
            i -= 1
            if i == 0:
                ret.add(BPD.from_bruhat_path(tuple(bpath + [w3 * Permutation.w0(len(w) - i).shiftup(i)])))
                continue
            da_list = pull_out_var(i, w3)
            #bpath.append(w)
            for _, perm in da_list:
                stack.append((perm, bpath + [w3 * Permutation.w0(len(w) - i).shiftup(i)], i))
        return ret
            

    w = Permutation([1,3,4,2])
    n = 3
    ws = Permutation.all_permutations(n)
    #for w in ws:
    bpds = BPD.all_bpds(w, n)
    rcs = {BPD.from_rc_graph(rc) for rc in RCGraph.all_rc_graphs(w, n)}
    for rc in rcs:
        if rc not in bpds:
            print("Boopy")
            print("rc:")
            pretty_print(rc)
    for bpd in bpds:
        print("Paint")
        pretty_print(bpd)
            # print("bpds:")
            # for bpd in bpds:
            #     pretty_print(bpd)
            # sys.exit(1)
    # # bpd = BPD.rothe_bpd(w, len(w))
    # print("FAT")
    # # pretty_print(bpd)

    # for rc in RCGraph.all_rc_graphs(w, len(w.trimcode)):
    #     bpd = BPD.from_rc_graph(rc)
    #     for i in range(len(rc)):
    #         print("dababa")
    #         pretty_print(bpd)
    #         print
    #     print()

    sys.exit(0)

    n = 5
    m = 4
    
        

    # row_ws = [w for w in ws if len(w) <= m]
    # grass_ws = [w for w in ws if w.descents() == {m - 1}]
    # for w in row_ws:
    #     for grass in grass_ws:
    #         for rc in RCGraph.all_rc_graphs(w, m):
    #             for grass_rc in RCGraph.all_rc_graphs(grass, m):
    #                 the_squash = bpd_squash(rc, grass_rc)
                    
    #                 print("rc:")
    #                 pretty_print(rc)
    #                 print("grass_rc:")
    #                 pretty_print(grass_rc)
    #                 print("squash:")
    #                 pretty_print(the_squash)
    #                 the_squash2 = rc.squash_product(grass_rc)
    #                 print("squash2:")
    #                 pretty_print(the_squash2)
    #                 if the_squash != the_squash2:
    #                     print("Boopy")
    #                     assert the_squash.length_vector == the_squash2.length_vector
                    
     