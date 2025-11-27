from schubmult.localization.leibniz_extension import compute_product
from schubmult.products.extended_products import new_product_method
from schubmult.crystal.non_principal_actions import NonPrincipalRCGraph
from schubmult.analysis.useful_pairs import find_useful_pairs

def simple_extension_example():
    # Example of using the new product method with a dominant RC-graph
    dominant_rc = ...  # Obtain or define a dominant RC-graph
    non_dominant_rc = ...  # Obtain or define a non-dominant RC-graph

    # Compute the product using the new method
    product_result = new_product_method(dominant_rc, non_dominant_rc)
    print("Product of dominant and non-dominant RC-graphs:")
    print(product_result)

def non_principal_example():
    # Example of using the new product method with non-principal RC-graphs
    non_principal_rc1 = NonPrincipalRCGraph(...)  # Define a non-principal RC-graph
    non_principal_rc2 = NonPrincipalRCGraph(...)  # Define another non-principal RC-graph

    # Compute the product of two non-principal RC-graphs
    product_result = new_product_method(non_principal_rc1, non_principal_rc2)
    print("Product of two non-principal RC-graphs:")
    print(product_result)

def main():
    print("Running simple extension examples...")
    simple_extension_example()
    non_principal_example()

if __name__ == "__main__":
    main()