from schubmult.localization.leibniz_extension import compute_leibniz_extension
from schubmult.products.extended_products import compute_extended_product
from schubmult.crystal.non_principal_actions import NonPrincipalRCGraph
import unittest

class TestLeibnizExtension(unittest.TestCase):

    def setUp(self):
        # Setup code to initialize any necessary variables or objects
        self.rc_graph_a = NonPrincipalRCGraph(...)  # Replace with actual initialization
        self.rc_graph_dom = NonPrincipalRCGraph(...)  # Replace with actual initialization
        self.index = 1  # Example index for testing

    def test_leibniz_extension(self):
        # Test the Leibniz extension method
        result = compute_leibniz_extension(self.rc_graph_a, self.rc_graph_dom, self.index)
        expected_result = ...  # Define the expected result based on your logic
        self.assertEqual(result, expected_result)

    def test_extended_product(self):
        # Test the computation of extended products
        result = compute_extended_product(self.rc_graph_a, self.rc_graph_dom)
        expected_result = ...  # Define the expected result based on your logic
        self.assertEqual(result, expected_result)

    def test_invalid_inputs(self):
        # Test handling of invalid inputs
        with self.assertRaises(ValueError):
            compute_leibniz_extension(None, self.rc_graph_dom, self.index)
        with self.assertRaises(ValueError):
            compute_extended_product(self.rc_graph_a, None)

if __name__ == '__main__':
    unittest.main()