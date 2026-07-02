from schubmult import *
from schubmult.rings.combinatorial.forest_rc_ring import ForestRCGraphRing
from schubmult.rings.polynomial_algebra import *
from schubmult.combinatorics.indexed_forests import weak_composition_to_indfor

f = ForestRCGraphRing()
ff = f @ f


def _coerce_to_tensor_elem(x):
    """Coerce scalars / compatible ring elements into f@f elements."""
    if hasattr(x, "ring") and x.ring == ff:
        return x
    if hasattr(x, "ring") and hasattr(x, "items"):
        try:
            return ff.from_comp_ring(x)
        except Exception:
            pass
    # Treat plain scalars/sympy atoms as coefficients on tensor identity.
    return ff.from_dict({ff.zero_monom: x})

def leibniz_forest(rc1, rc2, index):    
    return f(rc1).forest_trim(index) @ f(rc2).quasi_shift(index + 1) + f(rc1).quasi_shift(index) @ f(rc2).forest_trim(index)
    
def forest_trim_tensor(tensor_elem, index_list):
    tensor_elem = _coerce_to_tensor_elem(tensor_elem)
    if len(index_list) == 0:
        return tensor_elem
    index = index_list[-1]
    rest = index_list[:-1]
    result = ff.zero
    for (rc1, rc2), coeff in tensor_elem.items():
        result += coeff * leibniz_forest(rc1, rc2, index)
    return forest_trim_tensor(result, rest)
        

def comp_to_word(comp):
    if len(comp) == 0:
        return []
    word = []
    forest = weak_composition_to_indfor(comp)
    while forest.trim_descents:
        desc = forest.trim_descents[0]
        word = [desc, *word]
        forest = forest.trim_descent(desc)
    return word

if __name__ == "__main__":
    import sys
    import itertools

    n = int(sys.argv[1])
    comps = [perm.pad_code(n -  1) for perm in Permutation.all_permutations(n)]
    for comp1, comp2 in itertools.product(comps, repeat=2):
        rc1 = f.forest_poly(comp1)
        rc2 = f.forest_poly(comp2)
        prd = Forest(*comp1) * Forest(*comp2)
        pork = rc1 @ rc2
        for end_comp, val in prd.items():
            varnish = forest_trim_tensor(pork, comp_to_word(end_comp))
            #assert len(varnish) == 1, f"Error: Leibniz rule produced multiple terms for {comp1=} {comp2=} {end_comp=} {val=} {varnish=}"
            assert (sum(varnish.values())-val)==0, f"Error: Leibniz rule failed for {comp1=} {comp2=} {end_comp=} {val=} {varnish=}"
        print("Socks")