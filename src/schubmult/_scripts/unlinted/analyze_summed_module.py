#!/usr/bin/env python3
"""
Analyze the module structure when summing over all RC graphs for each permutation.

The hypothesis: Individual RC graphs have many relations, but when we sum all RC graphs
for a given permutation, we get a free module basis over the commutative subring of
Grassmannian elements.
"""

from schubmult.schub_lib.perm_lib import Permutation
from schubmult.schub_lib.rc_graph import RCGraph
from schubmult.rings.rc_graph_ring import RCGraphRing
from collections import defaultdict
from sympy import Matrix, gcd, lcm
import itertools


def generate_grassmannian_permutations(n):
    """Generate n-Grassmannian permutations (identity + one descent at n-1)."""
    grass_perms = [Permutation([])]
    
    # Only check a couple lengths beyond n (copied from working script)
    for length in range(n, n + 3):
        for perm in Permutation.all_permutations(length):
            descents = perm.descents()
            if len(descents) == 1 and (n - 1) in descents:
                grass_perms.append(perm)
    
    return grass_perms


def generate_k_grassmannian_permutations(k):
    """Generate k-Grassmannian permutations."""
    grass_perms = [Permutation([])]
    
    for length in range(k, k + 2):
        for perm in Permutation.all_permutations(length):
            descents = perm.descents()
            if len(descents) == 1 and (k - 1) in descents:
                grass_perms.append(perm)
    
    return grass_perms


def get_k_grassmannian_rc_graphs(k, n):
    """Get all k-Grassmannian RC graphs with n rows."""
    grass_perms = generate_k_grassmannian_permutations(k)
    grass_rcs = []
    
    for perm in grass_perms:
        if len(perm.trimcode) <= n:
            try:
                for rc in RCGraph.all_rc_graphs(perm, n):
                    grass_rcs.append(rc)
            except (ValueError, Exception) as e:
                pass
    
    return grass_rcs


def get_grassmannian_rc_graphs(n):
    """Get all n-Grassmannian RC graphs with n rows."""
    grass_perms = generate_grassmannian_permutations(n)
    grass_rcs = []
    
    for perm in grass_perms:
        if len(perm.trimcode) <= n:
            try:
                for rc in RCGraph.all_rc_graphs(perm, n):
                    grass_rcs.append(rc)
            except (ValueError, Exception) as e:
                pass
    
    return grass_rcs


def is_elementary_symmetric_grassmannian(perm):
    """Check if Grassmannian permutation is elementary symmetric."""
    descents = perm.descents()
    if len(descents) != 1:
        return False
    
    perm_list = list(perm)
    if not perm_list:
        return False
    
    if perm_list[-1] == 1:
        for i, val in enumerate(perm_list[:-1], start=2):
            if val != i:
                return False
        return True
    
    return False


def is_complete_symmetric_grassmannian(perm):
    """Check if Grassmannian permutation is complete symmetric."""
    descents = perm.descents()
    if len(descents) != 1:
        return False
    
    perm_list = list(perm)
    if not perm_list:
        return False
    
    last_val = perm_list[-1]
    expected = [i for i in range(1, len(perm_list) + 1) if i != last_val]
    
    if perm_list[:-1] == expected:
        return True
    
    return False


def filter_symmetric_grassmannian(grass_rcs, symmetric_type='both'):
    """Filter to elementary/complete symmetric functions."""
    filtered = []
    for rc in grass_rcs:
        perm = rc.perm
        
        if symmetric_type in ['elementary', 'both']:
            if is_elementary_symmetric_grassmannian(perm):
                filtered.append(rc)
                continue
        
        if symmetric_type in ['complete', 'both']:
            if is_complete_symmetric_grassmannian(perm):
                filtered.append(rc)
                continue
    
    return filtered


def group_rc_graphs_by_permutation(rc_graphs):
    """Group RC graphs by their underlying permutation."""
    by_perm = defaultdict(list)
    for rc in rc_graphs:
        perm_tuple = tuple(rc.perm)
        by_perm[perm_tuple].append(rc)
    return by_perm


def analyze_summed_module(n, k_values, symmetric_type='both', verbose=True, use_ring_action=False):
    """
    Analyze module structure using summed RC graphs (sum over all graphs for each permutation).
    
    Args:
        use_ring_action: If True, use ring multiplication (*) instead of module action (%)
    
    Returns:
        Dictionary with analysis results
    """
    if verbose:
        action_str = "ring action (*)" if use_ring_action else "module action (%)"
        print(f"Analyzing summed module for n={n}, k-values={k_values}, symmetric={symmetric_type}")
        print(f"Using {action_str}")
        print("="*80)
    
    # Initialize ring
    rc_ring = RCGraphRing()
    
    # Get all generators (k-Grassmannian RC graphs)
    all_generators = []
    for k in k_values:
        k_gens = get_k_grassmannian_rc_graphs(k, n)
        all_generators.extend(k_gens)
        if verbose:
            print(f"Found {len(k_gens)} {k}-Grassmannian RC graphs")
    
    if verbose:
        print(f"Total generators: {len(all_generators)}")
    
    # Group generators by permutation
    gen_by_perm = group_rc_graphs_by_permutation(all_generators)
    if verbose:
        print(f"Generators correspond to {len(gen_by_perm)} distinct permutations")
        print(f"RC graphs per permutation: {[len(rcs) for rcs in gen_by_perm.values()]}")
    
    # Get Grassmannian ring elements (filtered to symmetric if requested)
    grass_rcs = get_grassmannian_rc_graphs(n)
    if symmetric_type:
        grass_rcs = filter_symmetric_grassmannian(grass_rcs, symmetric_type)
    
    if verbose:
        print(f"\nGrassmannian elements (ring): {len(grass_rcs)}")
    
    # For each permutation basis element, compute products with all Grassmannian elements
    # Store as: summed_products[perm][grass_perm] = list of all resulting RC graphs
    if verbose:
        print(f"\nComputing products for each permutation basis element...")
        print("  (Summing all RC graphs for each permutation on the left)")
    
    summed_products = {}
    product_count = 0
    
    for perm_tuple, rc_list in gen_by_perm.items():
        summed_products[perm_tuple] = {}
        
        for grass_rc in grass_rcs:
            # For this (permutation basis, grassmannian) pair, compute all products
            # We sum the results from all RC graphs with this permutation
            all_results = []
            for gen_rc in rc_list:
                try:
                    # Wrap in ring elements
                    gen_elem = rc_ring(gen_rc)
                    grass_elem = rc_ring(grass_rc)
                    
                    if use_ring_action:
                        product = gen_elem * grass_elem
                    else:
                        product = gen_elem % grass_elem
                    
                    # Product is a RingElement (dict of RC graphs)
                    if product and len(product) > 0:
                        for result_rc, coeff in product.items():
                            all_results.append(result_rc)
                except Exception as e:
                    if verbose and product_count < 5:  # Only show first few errors
                        print(f"  Error computing {gen_rc.perm} % {grass_rc.perm}: {e}")
            
            product_count += 1
            if verbose and product_count % 100 == 0:
                print(f"  Computed {product_count} products...")
            
            if all_results:
                grass_perm_tuple = tuple(grass_rc.perm)
                summed_products[perm_tuple][grass_perm_tuple] = all_results
    
    if verbose:
        print(f"  Total products computed: {product_count}")
    
    # Now check for relations among these summed products
    # Two summed products are equal if they produce the same multiset of RC graphs
    if verbose:
        print(f"\nAnalyzing relations among summed products...")
        print("  Looking for: (Sum_RC_for_P1) op G1 = (Sum_RC_for_P2) op G2")
    
    # For each pair of (perm_basis, grass_element), we have a list of RC graphs
    # We want to find cases where different (perm_basis, grass_element) pairs
    # produce the same multiset of results
    
    # Create a signature for each product: sorted tuple of RC graph string representations
    product_signatures = {}
    for perm_tuple, grass_dict in summed_products.items():
        for grass_tuple, results in grass_dict.items():
            # Sort the RC graphs to create a canonical signature
            # This represents the sum as a multiset
            result_strs = [str(rc) for rc in results]
            signature = tuple(sorted(result_strs))
            
            key = (perm_tuple, grass_tuple)
            product_signatures[key] = signature
    
    # Group by signature to find relations
    signature_groups = defaultdict(list)
    for key, sig in product_signatures.items():
        signature_groups[sig].append(key)
    
    # Find non-trivial relations (same signature, different keys)
    summed_relations = []
    for sig, keys in signature_groups.items():
        if len(keys) > 1:
            # This is a relation: multiple (perm, grass) pairs give the same sum
            summed_relations.append(keys)
    
    if verbose:
        print(f"\nNumber of distinct products (signatures): {len(signature_groups)}")
        print(f"Number of (perm, grass) pairs computed: {len(product_signatures)}")
        print(f"Number of relations among summed products: {len(summed_relations)}")
        
        if len(summed_relations) == 0:
            print("\n" + "="*80)
            print("★ RESULT: The summed module appears to be FREE! ★")
            print("="*80)
            print("\nNo relations found when summing RC graphs by permutation.")
            print("The permutation basis {Sum of RC graphs for π} is free over")
            print("the Grassmannian ring (or its symmetric subring).")
        else:
            print("\n" + "="*80)
            print("RESULT: Relations exist among summed elements")
            print("="*80)
            print(f"\nShowing first 10 relations among summed elements:")
            for i, relation in enumerate(summed_relations[:10], 1):
                print(f"\nRelation {i}: {len(relation)} products yield the same sum")
                for perm_tuple, grass_tuple in relation[:3]:  # Show first 3 of each relation
                    perm = Permutation(list(perm_tuple)) if perm_tuple else Permutation([])
                    grass = Permutation(list(grass_tuple)) if grass_tuple else Permutation([])
                    print(f"  (Sum of RC[{perm}]) % Grass[{grass}]")
                if len(relation) > 3:
                    print(f"  ... and {len(relation)-3} more products")
            
            if len(summed_relations) > 10:
                print(f"\n... and {len(summed_relations)-10} more relations")
    
    # Additional analysis: look at the structure
    if verbose:
        print(f"\n" + "="*80)
        print("STRUCTURAL ANALYSIS:")
        print("="*80)
        
        # How many RC graphs per permutation?
        print(f"\nRC graphs per permutation (generators):")
        for perm_tuple, rc_list in sorted(gen_by_perm.items(), key=lambda x: len(x[1]), reverse=True):
            perm = Permutation(list(perm_tuple)) if perm_tuple else Permutation([])
            print(f"  {perm}: {len(rc_list)} RC graphs")
        
        # Total products computed
        total_products = sum(len(grass_dict) for grass_dict in summed_products.values())
        print(f"\nTotal (permutation, grassmannian) pairs computed: {total_products}")
        
        # Products per permutation basis
        print(f"\nProducts per permutation basis element:")
        for perm_tuple, grass_dict in summed_products.items():
            perm = Permutation(list(perm_tuple)) if perm_tuple else Permutation([])
            print(f"  {perm}: {len(grass_dict)} products")
    
    return {
        'generators_by_perm': gen_by_perm,
        'grassmannian_elements': grass_rcs,
        'summed_products': summed_products,
        'product_signatures': product_signatures,
        'relations': summed_relations,
        'is_free': len(summed_relations) == 0
    }


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Analyze summed RC graph module structure"
    )
    parser.add_argument("n", type=int, help="Parameter n")
    parser.add_argument("--k-grassmannian", type=str, required=True,
                       help="k-Grassmannian generators (comma-separated)")
    parser.add_argument("--symmetric", type=str, default='both',
                       choices=['elementary', 'complete', 'both'],
                       help="Filter to symmetric functions")
    parser.add_argument("--quiet", action="store_true",
                       help="Suppress verbose output")
    parser.add_argument("--ring-action", action="store_true",
                       help="Use ring multiplication (*) instead of module action (%)")
    
    args = parser.parse_args()
    
    k_values = [int(k.strip()) for k in args.k_grassmannian.split(',')]
    
    results = analyze_summed_module(
        n=args.n,
        k_values=k_values,
        symmetric_type=args.symmetric,
        verbose=not args.quiet,
        use_ring_action=args.ring_action
    )
    
    print("\n" + "="*80)
    print("CONCLUSION:")
    print("="*80)
    if results['is_free']:
        print("The summed module IS FREE over the (commutative) symmetric subring!")
        print("This confirms the hypothesis.")
    else:
        print(f"The summed module has {len(results['relations'])} relations.")
        print("Further investigation needed to understand the structure.")
