from schubmult import *
from sympy import pretty_print
if __name__ == "__main__":
    my_list_rc_2 = [RootTableau.from_rc_graph(rc) for rc in RCGraph.all_rc_graphs(Permutation.ref_product(4,5,1,4,6,2,3,2,5),6)]
    rt = my_list_rc_2[0]
    print("Original tableau:")
    pretty_print(rt)
    print(rt.reduced_word)
    delbox = (1,0)
    print(f"After deleting box {delbox}:")
    drt = rt.delete_box(*delbox)
    pretty_print(drt)
    print(drt.reduced_word)