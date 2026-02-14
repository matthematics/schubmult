from schubmult import *

P = PlacticAlgebra()
N = NilPlacticAlgebra()

def cauchy_kernel(index, max_size):
    elem = (P @ N).one
    for i in range(index):
        elem *= (P @ N).one + (P(Plactic(((index - i,),)))) @ N(NilPlactic(((max_size - i,),)))
    return elem
    
if __name__ == "__main__":
    import sys

    n = int(sys.argv[1])
    #perms = Permutation.all_permutations(n)
    full_kernel = (P @ N).one
    for i in range(n):
        full_kernel =  full_kernel * cauchy_kernel(i, n - 1) 
    dct = {}
    for (p, an), coeff in full_kernel.items():
        dct[an] = dct.get(an, P.zero) + coeff*P(p)

    for nil, coeff in dct.items():
        print(f"{nil.row_word}:\n{coeff}")
    