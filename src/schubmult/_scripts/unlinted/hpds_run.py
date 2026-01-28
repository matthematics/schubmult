from schubmult import *

if __name__ == "__main__":
    rc = RCGraph.random_rc_graph(uncode([5,6,0,0,2,3,5]))
    print(HPD.from_rc_graph(rc))
    print(HPD.from_bpd(BPD.from_rc_graph(rc)))