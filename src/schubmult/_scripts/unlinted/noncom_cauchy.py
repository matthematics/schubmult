from schubmult import *
from schubmult.utils.perm_utils import artin_sequences
from schubmult.rings.polynomial_algebra import *
from sympy import pretty_print

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)

    stairs = list(range(n - 1, 0, -1))
    g = GrassTensorAlgebra()
    r = RCGraphRing()
    bacon_elem = 0
    for seq in artin_sequences(n - 1):
        elem = g.one
        for i, a in enumerate(seq):
            elem *= g.elem_sym(stairs[i] - a, n - 1 - i)
        #the_new_elem = elem.to_rc_graph_ring_element().resize(n - 1)
        
        for tensor, coeff in elem.items():
            bacon_elem += coeff * g(tensor) @ g(tensor).to_rc_graph_ring_element().resize(n - 1) @ PA(*seq).change_basis(SchubertPolyBasis)
    pretty_print(bacon_elem)