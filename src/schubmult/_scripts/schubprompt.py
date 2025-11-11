import os

from prompt_toolkit import prompt
from prompt_toolkit.history import FileHistory
from prompt_toolkit.shortcuts import PromptSession
from sympy import pretty_print

from schubmult import *

# Define a history file for storing past commands
history_file = os.path.expanduser('~/.schubert_shell_history')

# Create a PromptSession with file history

session = PromptSession(history=FileHistory(history_file))


def main():
    my_list_rc_2 = [RootTableau.from_rc_graph(rc) for rc in RCGraph.all_rc_graphs(Permutation.ref_product(4,5,1,4,6,2,3,2,5),6)]
    rt = my_list_rc_2[0]
    print("Original tableau:")
    pretty_print(rt)
    print(rt.reduced_word)
    assert Permutation.ref_product(*rt.reduced_word) == Permutation.ref_product(4,5,1,4,6,2,3,2,5)

    while True:
            # Use the session's prompt method for input

        try:
            command = session.prompt('> ')
            if command.lower() == 'exit':
                break

            cmd, params = command.split("(", 1)
            params = params[:-1]
            if cmd == "ujdt":
                a, b = params.split(",")
                a, b = int(a), int(b)
                rt = rt.up_jdt_slide(a, b)
            elif cmd == "djdt":
                a, b = params.split(",")
                a, b = int(a), int(b)
                rt = rt.down_jdt_slide(a, b)
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
            elif cmd == "eval":
                eval(params)
            else:
                print(f"Invalid command: {cmd}")
        except EOFError:  # Ctrl-D
            break
        except KeyboardInterrupt:  # Ctrl-C
            print("Operation cancelled.")
            continue
        except Exception:
            import traceback
            print("Error occurred")
            traceback.print_exc()


        pretty_print(rt)
        print(rt.reduced_word)



if __name__ == '__main__':
    main()

