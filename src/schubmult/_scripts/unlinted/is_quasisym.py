from schubmult import *
from schubmult.symbolic import *
from itertools import combinations

def get_composition_and_indices(exponents):
    """
    Extracts the non-zero exponents (composition) and their positions (indices).
    Example: (2, 0, 1) -> Composition: (2, 1), Indices: (0, 2)
    """
    composition = []
    indices = []
    for i, exp in enumerate(exponents):
        if exp > 0:
            composition.append(exp)
            indices.append(i)
    return tuple(composition), tuple(indices)

def is_quasisymmetric(monomials, n):
    """
    Checks if a polynomial is quasisymmetric in n variables.
    
    :param monomials: A list of tuples (coefficient, exponent_tuple)
                      e.g., [(1, (2, 1, 0)), (1, (2, 0, 1))]
    :param n: The number of variables (e.g., 3 for x1, x2, x3)
    :return: True if quasisymmetric, False otherwise.
    """
    # map: composition -> { index_tuple: coefficient }
    comp_map = {}
    
    for exps, coeff in monomials.items():
        if coeff == 0: continue
        
        # Ensure the exponent tuple matches the number of variables
        if len(exps) < n:
            exps = exps + (0,) * (n - len(exps))
        elif len(exps) > n:
            # If the polynomial has terms with variables > n, it's not in n variables
            if any(e > 0 for e in exps[n:]):
                return False
            exps = exps[:n]

        comp, idxs = get_composition_and_indices(exps)
        
        if comp not in comp_map:
            comp_map[comp] = {}
        comp_map[comp][idxs] = coeff
        
    for comp, found_indices in comp_map.items():
        k = len(comp)
        # Generate all strictly increasing index sequences of length k from {0, ..., n-1}
        all_possible_indices = list(combinations(range(n), k))
        
        # 1. Check if every required index sequence exists in the polynomial
        if len(found_indices) != len(all_possible_indices):
            return False
            
        # 2. Check if all coefficients for this composition are identical
        coeffs = list(found_indices.values())
        first_coeff = coeffs[0]
        if not all(c == first_coeff for c in coeffs):
            return False
            
    return True

def monomial_quasisym(comp, length, genset):
    if any(c == 0 for c in comp):
        return S.Zero
    if len(comp) == length:
        return expand_seq(comp, genset)
    if len(comp) == 0:
        return S.One
    ret = S.Zero
    
    ret += monomial_quasisym(comp, length - 1, genset) + monomial_quasisym(comp[:-1], length - 1, genset) * genset[length]**comp[-1]
    return ret

if __name__ == "__main__":
    import sys
    from schubmult.rings.variables import genset_dict_from_expr
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)

    for perm in perms:
        num_vars = len(perm.trimcode)
        poly = Sx(perm).expand()
        dct = genset_dict_from_expr(poly, Sx.genset)
        print(poly)
        print(f"{perm.trimcode}: {is_quasisymmetric(dct, num_vars)}")
        spittoon = monomial_quasisym(perm.trimcode, num_vars+5, Sx.genset)
        print(f"Spittoon: {spittoon}")
        print(is_quasisymmetric(genset_dict_from_expr(spittoon, Sx.genset), num_vars+5))