from schubmult import *
import numpy as np

if __name__ == "__main__":
    rc = RCGraph.random_rc_graph(uncode([5,6,0,0,2,3,5]))

    rc1, rc2 = rc.rowrange(0, 4), rc.rowrange(4)
    bpd2 = BPD.from_rc_graph(rc2)
    
    #bpd2 = bpd2.delete_top_row()
    print("RC Graph:")
    print(rc)
    print("RC1:")
    print(rc1)
    print("BPD2:")
    print(bpd2)
    #assert rc1.perm * (bpd2.perm.shiftup(len(rc1))) == rc.perm
    combined = HPD.concat(rc1, bpd2)
    print("Combined HPD:")
    print(combined)
    print(f"{~(rc.perm)}")
    combined = combined.toggle_bottom_row()
    print("After swap:")
    print(combined)
    # # rc_bottom, bpd_top = rc.vertical_cut(4)
    # # bpd_top = BPD.from_rc_graph(bpd_top)
    # # print("RC Graph:")
    # # print(rc)
    # # print("RC Bottom:")
    # # print(rc_bottom)
    # # print("BPD Top:")
    # # print(bpd_top)
    # # combined = HPD.concat(bpd_top, rc_bottom)
    # # print("Combined HPD:")
    # # print(combined)
    # # rc = rc.resize(len(rc.perm))
    # # RCB = HPD.from_bpd(BPD.from_rc_graph(rc))
    # # for row in range(RCB.rows):
    # #     print(RCB[row, :])
    # #     for col in range(RCB.cols):
    # #         print(f"{(row,col)} {RCB.pipe_source_labels(row, col)=} {RCB[row, col]}")
    # # print(RCB)
    # # #print(RCB.length_vector)
    # # BBB = HPD.from_bpd(BPD.from_rc_graph(rc))
    # # print(BBB)
    # # new_grid = np.concatenate((RCB._grid[:4, :], BBB._grid[4:, :]), axis=0)
    # # print(HPD(new_grid, (0,0,0,0,1,1,1)))
    # # print(BBB)
    # # print(BBB.length_vector)
    # # print("THE TOG")
    # # print(BBB.toggle_bottom_row())
    # # grid = np.array([[HPDTile.CROSS, HPDTile.BUMP, HPDTile.CROSS, HPDTile.ELBOW_NW],
    # #                  [HPDTile.VERT, HPDTile.ELBOW_NE, HPDTile.CROSS, HPDTile.HORIZ],
    # #                  [HPDTile.ELBOW_NW, HPDTile.ELBOW_SE, HPDTile.ELBOW_NW, HPDTile.BLANK],
    # #                  [HPDTile.HORIZ, HPDTile.ELBOW_NW, HPDTile.BLANK, HPDTile.BLANK]])
    # # hpd = HPD(grid, (0,1,0,0))
    # # print(hpd)
    
    # # # grid = np.array([[HPDTile.CROSS, HPDTile.BUMP, HPDTile.BUMP, HPDTile.ELBOW_NW],
    # # #                  [HPDTile.VERT, HPDTile.ELBOW_NE, HPDTile.CROSS, HPDTile.HORIZ],
    # # #                  [HPDTile.ELBOW_NW, HPDTile.ELBOW_SE, HPDTile.ELBOW_NW, HPDTile.BLANK],
    # # #                  [HPDTile.BLANK, HPDTile.ELBOW_NE, HPDTile.HORIZ, HPDTile.HORIZ]])
    # # grid = np.array([[HPDTile.BUMP, HPDTile.BUMP, HPDTile.ELBOW_NW],
    # #                  [HPDTile.CROSS, HPDTile.ELBOW_NW, HPDTile.BLANK],
    # #                  [HPDTile.ELBOW_NE, HPDTile.HORIZ, HPDTile.HORIZ]])
    # # hpd = HPD(grid, (0,0,1))
    # hpd = HPD.from_rc_graph(rc)
    # print(hpd)
    # hpd = hpd.toggle_bottom_row()
    # print(hpd)
    # print("swapping row 11 and 12")
    # hpd = hpd.swap_rows(11)
    # print(hpd)
    # hpd = hpd.toggle_bottom_row()
    # print(hpd)
    # print("swapping row 10 and 11")
    # hpd = hpd.swap_rows(10)
    # print(hpd)
    # print("swapping row 11 and 12")
    # hpd = hpd.swap_rows(11)
    # print(hpd)
    # hpd = hpd.toggle_bottom_row()
    # print(hpd)
    # hpd = hpd.swap_rows(9).swap_rows(10).swap_rows(11)
    # print(hpd)
    # hpd = hpd.toggle_bottom_row()
    # print(hpd)
    # hpd = hpd.swap_rows(8).swap_rows(9).swap_rows(10).swap_rows(11)
    # print(hpd)
    # hpd = hpd.toggle_bottom_row()
    # print(hpd)
    # hpd = hpd.swap_rows(7).swap_rows(8).swap_rows(9).swap_rows(10).swap_rows(11)
    # print(hpd)
    # hpd = hpd.toggle_bottom_row()
    # print(hpd)
    # hpd = hpd.swap_rows(6).swap_rows(7).swap_rows(8).swap_rows(9).swap_rows(10).swap_rows(11)
    # print(hpd)
    # hpd = hpd.toggle_bottom_row()
    # print(hpd)

    # bpd = hpd.to_bpd(6, 12)
    # # for row in range(hpd.rows):
    # #     for col in range(hpd.cols):
    # #         bungee = set([a for a in hpd.pipe_source_labels(row, col).values() if a is not None])
    # #         assert hpd[row, col] == HPDTile.BLANK or any(a<=len(hpd.perm.trimcode) for a in bungee) or not hpd.is_weighty_position(row, col)
    # print(bpd)
    # print("pants")