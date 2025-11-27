from schubmult.localization.leibniz_extension import compute_product
from schubmult.localization.weight_tracker import WeightTracker
from schubmult.products.extended_products import ExtendedProductCalculator
from schubmult.crystal.non_principal_actions import NonPrincipalActions
from schubmult.analysis.useful_pairs import find_useful_pairs

def verify_weights(rc_graphs, dominant_rcs):
    """
    Verify the weights of products computed from RC-graphs.
    
    Parameters:
    - rc_graphs: List of RC-graphs to compute products from.
    - dominant_rcs: List of dominant RC-graphs for comparison.
    
    Returns:
    - A dictionary with RC-graph permutations as keys and their verified weights as values.
    """
    weight_tracker = WeightTracker()
    verified_weights = {}

    for rc in rc_graphs:
        for dominant_rc in dominant_rcs:
            product = compute_product(rc, dominant_rc)
            weight_tracker.track_weight(product)
            verified_weights[rc.perm.trimcode] = product.crystal_weight

    return verified_weights

if __name__ == "__main__":
    # Example usage
    from schubmult import RCGraph, Permutation

    # Define some RC-graphs and dominant RC-graphs for testing
    rc_graphs = [RCGraph.principal_rc(Permutation([1, 2, 3]), 2)]
    dominant_rcs = [RCGraph.principal_rc(Permutation([1, 2, 3]), 2)]

    verified_weights = verify_weights(rc_graphs, dominant_rcs)
    for perm, weight in verified_weights.items():
        print(f"RC-graph {perm} has verified weight: {weight}")