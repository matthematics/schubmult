from schubmult import *
from schubmult.symbolic import prod
from schubmult.abc import *
from schubmult.utils.perm_utils import artin_sequences

if __name__ == "__main__":
    import sys

    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)
    w0 = Permutation.w0(n)
    dct = {}
    big_monom = Sx(w0).expand()
    for seq in artin_sequences(n - 1):
        the_monom = big_monom * DSx([]) * prod([H(seq[a], 1, [x[a+1]], tuple([yy**-1 for yy in y[1:n+5]])) for a in range(len(seq))] )
        for permo, coeff in the_monom.items():
            dct[permo] = dct.get(permo, 0) + coeff*FA(*seq)
    for permo, cc in dct.items():
        print(f"{permo}: {cc}")