from schubmult import RCGraph, RCGraphRing, Permutation

def compute_new_product(rc_a, rc_b, index):
    """
    Compute a new product of RC-graphs that is not limited to the form of
    (any RC) % (dominant RC). This method utilizes insights from divided
    differences and the Leibniz extension.

    Parameters:
    - rc_a: The first RC-graph (can be any RC).
    - rc_b: The second RC-graph (can be any RC).
    - index: The index for the divided difference operation.

    Returns:
    - A new RC-graph representing the product of rc_a and rc_b.
    """
    rc_ring = RCGraphRing()

    # Compute the product of the two RC-graphs
    product = rc_ring(rc_a) % rc_ring(rc_b)

    # Check the divided difference for the product
    divdiff_result = product.divdiff(index)

    if len(divdiff_result) == 0:
        # If the divided difference is zero, we can apply the Leibniz extension
        # to compute the product in a different way
        s_i_a = rc_ring(rc_a).divdiff(index)
        s_i_b = rc_ring(rc_b).divdiff(index)

        if len(s_i_a) > 0 and len(s_i_b) > 0:
            # If both divided differences are non-zero, we can combine them
            new_terms = []
            for term_a, coeff_a in s_i_a.items():
                for term_b, coeff_b in s_i_b.items():
                    new_terms.append((term_a % term_b, coeff_a * coeff_b))

            return new_terms
        else:
            # Handle cases where one of the divided differences is zero
            return []

    return divdiff_result


def analyze_useful_pairs(pairs):
    """
    Analyze useful pairs of RC-graphs for extending products.

    Parameters:
    - pairs: A list of tuples containing pairs of RC-graphs.

    Returns:
    - A list of results from the new product computations.
    """
    results = []
    for rc_a, rc_b, index in pairs:
        result = compute_new_product(rc_a, rc_b, index)
        results.append((rc_a.perm.trimcode, rc_b.perm.trimcode, index, result))
    
    return results


def find_useful_pairs(n):
    """
    Identify useful pairs of RC-graphs for extending products.

    Parameters:
    - n: The size of the permutations.

    Returns:
    - A list of useful pairs of RC-graphs.
    """
    perms = Permutation.all_permutations(n)
    useful_pairs = []

    for perm_a in perms:
        rc_a = RCGraph.principal_rc(perm_a, n-1)
        for perm_b in perms:
            rc_b = RCGraph.principal_rc(perm_b, n-1)
            for index in range(1, n):
                if rc_ring(rc_a).divdiff(index) == 0 and rc_ring(rc_b).divdiff(index) != 0:
                    useful_pairs.append((rc_a, rc_b, index))

    return useful_pairs


if __name__ == "__main__":
    n = 3  # Example size
    pairs = find_useful_pairs(n)
    results = analyze_useful_pairs(pairs)

    for rc_a, rc_b, index, result in results:
        print(f"Product of {rc_a} and {rc_b} at index {index}: {result}")