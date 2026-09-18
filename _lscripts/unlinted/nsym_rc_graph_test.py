from schubmult import *
from symengine import Symbol
from sympy import pretty_print

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


if __name__ == "__main__":
    import sys
    import itertools
    import gc
    from sympy import init_printing
    init_printing(wrap_line=False)
    n = int(sys.argv[1])
    
    # Generate connected permutations as a list (needed for combinations)
    all_perms = Permutation.all_permutations(n)
    connected_perms = [perm for perm in all_perms if 0 not in perm.trimcode]
    del all_perms  # Free memory immediately
    gc.collect()
    
    J = FreeAlgebra(JBasis)
    M = QSym()
    
    print(f"Computing composition dictionaries for {len(connected_perms)} permutations...")
    the_comp_dicts = {}
    
    # Compute comp_dicts one at a time
    for idx, perm in enumerate(connected_perms):
        comp_dict = {}
        for i in range(perm.inv - len(perm.trimcode) + 1):
            thisperm = perm.shiftup(i)
            # Process RC graphs one at a time instead of storing all
            for rc in RCGraph.all_rc_graphs(thisperm, len(thisperm.trimcode)):
                last_zero = max((idx_val for idx_val, val in enumerate(thisperm.trimcode) if val == 0), default=-1)
                if not all(val == 0 for val in thisperm.trimcode[:last_zero + 1]):
                    continue
                if 0 in rc.length_vector:
                    continue
                comp_dict[rc.length_vector] = comp_dict.get(rc.length_vector, 0) + 1
        
        # Only print if comp_dict is non-empty to reduce output
        if comp_dict:
            quas = sum([coeff * M(*comp) for comp, coeff in comp_dict.items()]) 
            print(f"[{idx+1}/{len(connected_perms)}] QSchub{perm.trimcode} = {str(quas)}")
        
        the_comp_dicts[perm] = comp_dict
        
        # Periodic garbage collection
        if (idx + 1) % 10 == 0:
            gc.collect()
    
    print(f"\nVerifying orthogonality for {len(connected_perms)*(len(connected_perms)+1)//2} pairs...")
    
    # Cache J basis conversions to avoid recomputing
    J_words_cache = {}
    for perm in connected_perms:
        J_words_cache[perm] = J(*perm.trimcode).change_basis(WordBasis)
    
    pair_count = 0
    total_pairs = len(connected_perms) * (len(connected_perms) + 1) // 2
    
    for perm1, perm2 in itertools.combinations_with_replacement(connected_perms, 2):
        pair_count += 1
        J_words1 = J_words_cache[perm1]
        J_words2 = J_words_cache[perm2]
        
        the_sum = 0
        words = set(J_words1.keys()) | set(J_words2.keys())
        for word in words:
            nsym_coeff = J_words1.get(word, 0)
            nsym_coeff2 = J_words2.get(word, 0)
            quasisym_coeff2 = the_comp_dicts[perm2].get(word, 0)
            quasisym_coeff = the_comp_dicts[perm1].get(word, 0)
            total_coeff = nsym_coeff * nsym_coeff2 * quasisym_coeff * quasisym_coeff2
            if total_coeff != 0:
                if total_coeff != 1:
                    print(f"[{pair_count}/{total_pairs}] Nontrivial coeff: {total_coeff}")
            the_sum += total_coeff
        
        if perm1 == perm2:
            assert the_sum == 1, f"Expected 1 for {perm1} * {perm2}, got {the_sum}"
        else:
            assert the_sum == 0, f"Expected 0 for {perm1} * {perm2}, got {the_sum}"
        
        # Progress indicator and periodic GC
        if pair_count % 50 == 0:
            print(f"Progress: {pair_count}/{total_pairs} pairs verified")
            gc.collect()
    
    print(f"\nAll {total_pairs} pairs verified successfully!")
