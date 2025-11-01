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
                rt0 = rt.raising_operator(index)
                if rt0 is None:
                    print("Raising operator returned None")
                else:
                    rt = rt0
            elif cmd == "l":
                index = int(params)
                rt0 = rt.lowering_operator(index)
                if rt0 is None:
                    print("Lowering operator returned None")
                else:
                    rt = rt0
            else:
                print(f"Invalid command: {cmd}")
        except Exception as e:
            import traceback
            print("Error occurred")
            traceback.print_exc()
            

        pretty_print(rt)
        print(rt.reduced_word)
        