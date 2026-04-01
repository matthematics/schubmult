from schubmult import *
from schubmult.utils.perm_utils import artin_sequences
from sympy import pretty_print
# EE = FreeAlgebra(ElementaryBasis)
# CSx = FreeAlgebra(CompositionSchubertBasis)
# Key = FreeAlgebra(KeyBasis)
# ForestDual = FreeAlgebra(ForestBasis)
# FA = FreeAlgebra(WordBasis)

r = RCGraphRing()
# quasi_yam ring

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    check = set()
    for perm in Permutation.all_permutations(n):
        for rc in RCGraph.all_rc_graphs(perm, n - 1):
            if not rc.is_quasi_yamanouchi:
                continue
            the_tuple = (rc.perm_word, rc.length_vector)
            if the_tuple in check:
                pretty_print(rc)
                raise ValueError(f"Duplicate RC graph data found for {the_tuple}.")
                
            check.add(the_tuple)
    print("No duplicates found in RC graph data.")
