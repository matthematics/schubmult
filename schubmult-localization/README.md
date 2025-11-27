# README for Schubert Multiplier Localization

## Overview

The Schubert Multiplier Localization project provides a framework for computing products of RC-graphs, extending beyond traditional methods that are limited to dominant RCs. By leveraging insights from divided differences and the Leibniz extension, this project aims to facilitate the exploration of new product structures in the context of Schubert calculus.

## Features

- **Leibniz Extension**: Implements methods for extending products of RC-graphs, allowing for more flexible combinations.
- **Weight Tracking**: Ensures that weight preservation is maintained during product computations.
- **Product Caching**: Optimizes performance by caching products of RC-graphs to avoid redundant calculations.
- **Extended Product Computation**: Provides functionality for computing new products of RC-graphs that are not restricted to dominant forms.
- **Consistency Checks**: Verifies the consistency of product computations to ensure correctness.

## Installation

To install the necessary dependencies, run:

```
pip install -r requirements.txt
```

## Usage

### Basic Example

To compute products of RC-graphs, you can use the provided methods in the `src/products/extended_products.py` file. Here’s a simple example:

```python
from src.products.extended_products import compute_product

# Example RC-graphs
rc_graph_a = ...
rc_graph_b = ...

# Compute the product
product = compute_product(rc_graph_a, rc_graph_b)
```

### Running Tests

To ensure everything is functioning correctly, you can run the unit tests included in the `tests` directory:

```
pytest tests/
```

## Examples

The `examples` directory contains scripts demonstrating various functionalities:

- `simple_extensions.py`: Shows how to use the new product methods for simple cases.
- `non_principal_products.py`: Illustrates the capabilities of the new methods with non-principal RC-graphs.
- `weight_verification.py`: Verifies that the computed products maintain the expected weights.

## Contributing

Contributions are welcome! Please feel free to submit a pull request or open an issue for any enhancements or bug fixes.

## License

This project is licensed under the MIT License. See the LICENSE file for more details.