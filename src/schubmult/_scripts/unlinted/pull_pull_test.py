from schubmult import *
from schubmult.symbolic import *
from schubmult.abc import *
def is_decomposable(w):
    for i in range(1, len(w) - 1):
        coset, w_J = w.coset_decomp(*list(range(1, i + 1)), *list(range(i + 2, len(w))))
        if coset.inv == 0 and set(w_J.code[: i + 1]) != {0} and set(w_J.code[i + 2 :]) != {0}:
            return True
    return False

if __name__ == "__main__":
    import sys
    from schubmult.utils.schub_lib import pull_out_var

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    for perm in perms:
        for i in perm.descents(zero_indexed=False):
            summ = S.Zero
            for rc in RCGraph.all_rc_graphs(perm, len(perm.trimcode)):
                if len(rc[i - 1]) != 0:
                    continue 
                pullout = rc.pull_out_row(i)
                # print(f"{perm=} {i=} pullout")
                # print(pullout)
                # print(rc.perm.inv - pullout.perm.inv)
                
                # for pw, rc0 in pullout:
                #     summ += (x[i]**len(pw))*rc0.polyvalue(x[:i] + x[i+1:])
                
                print(f"RC:\n{rc}")
                print(tuple(pullout))
                assert pullout.perm.inv == perm.inv
                summ += (x[i] ** (perm.inv - pullout.perm.inv)) * pullout.polyvalue(x[:i] + x[i+1:])
            assert expand(summ - Sx(perm).expand().subs(x[i], S.Zero)) == S.Zero, f"Error: pull out variable mismatch for permutation {perm} at row {i}:\nComputed sum:\n{summ}\nExpected:\n{Sx(perm).expand().subs(x[i], S.Zero)}\n{pullout}"
            print(f"Permutation: {perm}, length: {i}, verified.")
