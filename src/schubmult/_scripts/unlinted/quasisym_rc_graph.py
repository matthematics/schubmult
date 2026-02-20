from schubmult import *
from symengine import Symbol, S
from sympy import pretty_print
from schubmult.symbolic.poly.variables import genset_dict_from_expr
from schubmult.rings.quasisymmetric_functions import monomial_quasisym
from itertools import combinations

r = RCGraphRing()
t = Symbol('t')

def act_nsym_rc(alpha, R):
    if len(alpha) == 0:
        return r(R)
    if 0 in alpha:
        raise ValueError("Compositions must have positive parts.")
    actor = alpha[-1]
    result = r(RCGraph.one_row(actor)) * act_nsym_rc(alpha[:-1], R)
    pruned_result = r.zero
    for rc, nsym_coeff in result.items():
        if rc.perm.inv > 0 and rc.perm.trimcode[0] == 0:
            pruned_result += nsym_coeff * t * r(rc)
        else:
            pruned_result += nsym_coeff * r(rc)
    return pruned_result

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


if __name__ == "__main__":
    import sys
    import itertools
    import gc
    from sympy import init_printing
    init_printing(wrap_line=False)
    n = int(sys.argv[1])
    m = int(sys.argv[2]) if len(sys.argv) > 2 else n
    
    # Generate connected permutations as a list (needed for combinations)
    perms = Permutation.all_permutations(n)
    #connected_perms = [perm for perm in all_perms if 0 not in perm.trimcode]
    # del all_perms  # Free memory immediately
    # gc.collect()
    
    # print(f"Computing composition dictionaries for {len(connected_perms)} permutations...")
    the_comp_dict = {}
    
    # Compute comp_dicts one at a time
    for perm in perms:
        if perm.inv > 0 and perm.trimcode[0] == 0:
            continue
        for i in range(m - len(perm.trimcode) + 1):
            thisperm = perm.shiftup(i)
            # last_zero = max((idx_val for idx_val, val in enumerate(thisperm.trimcode) if val == 0), default=-1)
            # if not all(val == 0 for val in thisperm.trimcode[:last_zero + 1]):
            #     continue
            # Process RC graphs one at a time instead of storing all
            for rc in RCGraph.all_rc_graphs(thisperm, len(thisperm.trimcode)):
                if any(val == 0 for val in rc.length_vector):
                    continue
                the_comp_dict[perm] = the_comp_dict.get(perm, S.Zero) + monomial_quasisym(rc.length_vector, m, Sx.genset)
        
        # Only print if comp_dict is non-empty to reduce output
        # if comp_dict:
        #     quas = sum([coeff * M(*comp) for comp, coeff in comp_dict.items()]) 
        #     print(f"[{idx+1}/{len(connected_perms)}] QSchub{perm.trimcode} = {str(quas)}")
        
        # the_comp_dicts[perm] = comp_dict
        
        # # Periodic garbage collection
        # if (idx + 1) % 10 == 0:
        #     gc.collect()
        if perm not in the_comp_dict:
            print(f"No RC graphs found for {perm.trimcode}")
            continue
        assert is_quasisymmetric(genset_dict_from_expr(the_comp_dict[perm], Sx.genset), m), f"Failed quasisymmetry for {perm.trimcode} {the_comp_dict[perm]}"
        if 0 in perm.trimcode:
            print(f"{perm.trimcode}: nontrivial quasi")
    
    