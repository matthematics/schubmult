# schubmult

Littlewood-Richardson coefficients for Schubert polynomials (ordinary, double,
quantum, and quantum double), with optional parabolic subgroup support.

- **[Module reference](modules/README.md)** — one page per module, generated from source docstrings.
- [GitHub repository](https://github.com/matthematics/schubmult)
- [PyPI package](https://pypi.org/project/schubmult/)

## Getting started

```bash
pip install schubmult
```

```python
from schubmult import Sx

# Product of two Schubert polynomials, indexed by permutations
print(Sx([3, 1, 2]) * Sx([2, 1, 3]))
```

See the [module reference](modules/README.md) for the full API, starting with
[`schubmult.combinatorics.permutation`](modules/schubmult/combinatorics/permutation.md)
(the core `Permutation` object) and
[`schubmult.rings.schubert`](modules/schubmult/rings/schubert/index.md)
(the `Sx`/`DSx` ring interface).
