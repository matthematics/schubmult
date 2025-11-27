from schubmult.localization.leibniz_extension import compute_leibniz_extension
from schubmult.crystal.non_principal_actions import apply_non_principal_action
from schubmult.analysis.useful_pairs import find_useful_pairs

def compute_product(rc_a, rc_b):
    """
    Compute the product of two RC-graphs, rc_a and rc_b, utilizing insights from
    divided differences and the Leibniz extension method. This method allows for
    products that are not limited to the form of (any RC) % (dominant RC).

    Parameters:
    - rc_a: The first RC-graph.
    - rc_b: The second RC-graph.

    Returns:
    - A new RC-graph representing the product of rc_a and rc_b.
    """
    # Check if rc_b is dominant
    if rc_b.perm.is_dominant:
        return rc_a % rc_b  # Use standard product if rc_b is dominant

    # Compute the Leibniz extension for the product
    leibniz_result = compute_leibniz_extension(rc_a, rc_b)

    # Apply any necessary non-principal actions to the result
    extended_result = apply_non_principal_action(leibniz_result)

    return extended_result

def validate_product(rc_a, rc_b, product):
    """
    Validate the computed product of two RC-graphs to ensure consistency
    and correctness according to defined product rules.

    Parameters:
    - rc_a: The first RC-graph.
    - rc_b: The second RC-graph.
    - product: The computed product of rc_a and rc_b.

    Returns:
    - A boolean indicating whether the product is valid.
    """
    # Implement validation logic based on product rules
    # This could involve checking weights, descents, and other properties
    # For now, we will return True as a placeholder
    return True

def compute_all_products(rc_graphs):
    """
    Compute products for all pairs of RC-graphs in the provided list.

    Parameters:
    - rc_graphs: A list of RC-graphs.

    Returns:
    - A list of tuples containing (rc_a, rc_b, product).
    """
    products = []
    for i, rc_a in enumerate(rc_graphs):
        for rc_b in rc_graphs[i + 1:]:
            product = compute_product(rc_a, rc_b)
            if validate_product(rc_a, rc_b, product):
                products.append((rc_a, rc_b, product))
    return products

def main():
    """
    Main function to demonstrate the computation of products.
    This can be expanded or modified for testing purposes.
    """
    # Example usage
    rc_graphs = []  # Populate with actual RC-graphs
    products = compute_all_products(rc_graphs)
    for rc_a, rc_b, product in products:
        print(f"Product of {rc_a.perm.trimcode} and {rc_b.perm.trimcode} yields {product.perm.trimcode}")

if __name__ == "__main__":
    main()