from schubmult.localization.leibniz_extension import compute_product
from schubmult.products.extended_products import ExtendedProduct
from schubmult.products.consistency_check import check_consistency
import unittest

class TestConsistency(unittest.TestCase):

    def setUp(self):
        # Setup code to initialize necessary objects or states for tests
        self.rc_graph_a = ExtendedProduct(...)  # Replace with actual initialization
        self.rc_graph_b = ExtendedProduct(...)  # Replace with actual initialization

    def test_product_consistency(self):
        # Test that the product of two RC-graphs maintains consistency
        product = compute_product(self.rc_graph_a, self.rc_graph_b)
        self.assertTrue(check_consistency(product), "Product does not maintain consistency.")

    def test_weight_preservation(self):
        # Test that weight is preserved in the product computation
        weight_a = self.rc_graph_a.crystal_weight
        weight_b = self.rc_graph_b.crystal_weight
        product = compute_product(self.rc_graph_a, self.rc_graph_b)
        self.assertEqual(product.crystal_weight, weight_a + weight_b, "Weight preservation failed.")

    def test_invalid_product(self):
        # Test that an invalid product raises an appropriate exception
        with self.assertRaises(ValueError):
            compute_product(self.rc_graph_a, None)  # Example of invalid input

if __name__ == '__main__':
    unittest.main()