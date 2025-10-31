from sympy import pretty_print

from schubmult import *

if __name__ == "__main__":
    my_list_rc_2 = [RootTableau.from_rc_graph(rc) for rc in RCGraph.all_rc_graphs(Permutation.ref_product(4,5,1,4,6,2,3,2,5),6)]
    rt = my_list_rc_2[0]
    print("Original tableau:")
    pretty_print(rt)
    print(rt.reduced_word)
    assert Permutation.ref_product(*rt.reduced_word) == Permutation.ref_product(4,5,1,4,6,2,3,2,5)
    while True:
        delbox = input("Enter box to add (row,col), q to quit: ")
        if delbox.lower() == 'q':
            break
        try:
            cmd, params = delbox.split("(")
            params = params[:-1]
            if cmd == "ujdt":
                a, b = params.split(",")
                a, b = int(a), int(b)
                rt = rt.up_jdt_slide(a, b)
            elif cmd == "r":
                index = int(params)
                rt = rt.raising_operator(index)
            elif cmd == "l":
                index = int(params)
                rt = rt.lowering_operator(index)
        except Exception:
            print(f"Error occurred: {e}")

        pretty_print(rt)
        print(rt.reduced_word)
        