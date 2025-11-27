from schubmult.localization.leibniz_extension import compute_leibniz_extension
from schubmult.products.extended_products import compute_non_principal_product
from schubmult.crystal.non_principal_actions import NonPrincipalRCGraph
from schubmult.analysis.useful_pairs import find_useful_pairs

def showcase_non_principal_products(n):
    """
    Showcase examples of products involving non-principal RC-graphs.
    
    This function demonstrates how to compute products of non-principal RC-graphs
    using the extended product methods and the insights gained from divided differences.
    """
    # Find useful pairs of RC-graphs for extension
    useful_pairs = find_useful_pairs(n)
    
    print(f"Found {len(useful_pairs)} useful pairs for n={n}.")

    for rc_a, rc_b in useful_pairs:
        print(f"\nComputing product for RC-graphs: {rc_a.perm.trimcode} and {rc_b.perm.trimcode}")
        
        # Compute the product using the extended product method
        product = compute_non_principal_product(rc_a, rc_b)
        
        print(f"Product of {rc_a.perm.trimcode} and {rc_b.perm.trimcode}:")
        print(f"  Resulting weight: {product.crystal_weight}")
        print(f"  Resulting permutation: {product.perm.trimcode}")

if __name__ == "__main__":
    n = 4  # Example size, can be adjusted
    showcase_non_principal_products(n)