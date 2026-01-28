from schubmult import *
import numpy as np

if __name__ == "__main__":
    rc = RCGraph.random_rc_graph(uncode([5,6,0,0,2,3,5]))
    RCB = HPD.from_rc_graph(rc)
    print(RCB)
    print(RCB.length_vector)
    BBB = HPD.from_bpd(BPD.from_rc_graph(rc))
    print(BBB)
    print(BBB.length_vector)
    grid = np.array([[HPDTile.CROSS, HPDTile.BUMP, HPDTile.CROSS, HPDTile.ELBOW_NW],
                     [HPDTile.VERT, HPDTile.ELBOW_NE, HPDTile.CROSS, HPDTile.HORIZ],
                     [HPDTile.ELBOW_NW, HPDTile.ELBOW_SE, HPDTile.ELBOW_NW, HPDTile.BLANK],
                     [HPDTile.HORIZ, HPDTile.ELBOW_NW, HPDTile.BLANK, HPDTile.BLANK]])
    hpd = HPD(grid, (0,1,0,0))
    print(hpd)
    grid = np.array([[HPDTile.CROSS, HPDTile.BUMP, HPDTile.CROSS, HPDTile.ELBOW_NW],
                     [HPDTile.VERT, HPDTile.ELBOW_NE, HPDTile.CROSS, HPDTile.HORIZ],
                     [HPDTile.ELBOW_NW, HPDTile.ELBOW_SE, HPDTile.ELBOW_NW, HPDTile.BLANK],
                     [HPDTile.BLANK, HPDTile.ELBOW_NE, HPDTile.HORIZ, HPDTile.HORIZ]])
    hpd = HPD(grid, (0,1,0,1))
    print(hpd)
