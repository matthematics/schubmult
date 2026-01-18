# schubmult

## Program and package for rapid computation of Littlewood-Richardson coefficients of Schubert polynomials, compliant with sympy/symengine (and hence indirectly Sage)

The main purpose of this python package is for executing scripts to compute coefficients of products of 
various types of Schubert polynomials. Coproducts can also be computed, as well as substitution of 
commuting difference operators for quantum double Schubert polynomials. Quantum multiplication also has 
parabolic subgroup support for equivariant/triple, computed via the Peterson-Woodward comparison theorem and generalizations as 
per Huang/Li.

[Docs to be hosted on Wiki](https://github.com/matthematics/schubmult/wiki/schubmult-home)

## To install dev version

```
pip install git+https://github.com/matthematics/schubmult.git
```

## RCGraph and BPD Functionality

The package implements two main combinatorial models for Schubert calculus:

- **RCGraph (Reduced Compatible Graphs):** Encodes reduced words for permutations as graphs, supporting crystal operations and algebraic manipulations.
- **BPD (Bumpless Pipe Dreams):** Represents tilings of an $n \times n$ grid with local tile rules, providing an alternative model for Schubert polynomials.

### Bijections and Conversions

There is a canonical bijection (Gao and Huang, 2017) between RCGraphs and BPDs for a given permutation and grid size:
- `BPD.from_rc_graph(rc_graph)`: Converts an RCGraph to a BPD using the inversion data.
- `BPD.to_rc_graph()`: Converts a BPD back to its corresponding RCGraph by extracting the reduced compatible sequence from the tile configuration.
These conversions are invertible up to normalization and grid size.

### Operations

RCGraph and BPD objects support:
- Enumeration for a given permutation and grid size
- Crystal operators (raising/lowering) and combinatorial mutations (currently BPDs only through the bijection)
- Conversion to algebraic elements in the Schubert and nilHecke rings
- Visualization and pretty-printing

### Example Usage

```python
from schubmult.schub_lib.rc_graph import RCGraph
from schubmult.schub_lib.bpd import BPD
import random

rc = random.choice(list(RCGraph.all_rc_graphs(Permutaiton([5,1,6,2,4,3]), 4)))
bpd = BPD.from_rc_graph(rc)
rc2 = bpd.to_rc_graph()
assert rc2.perm == rc.perm
```

See the [Wiki](https://github.com/matthematics/schubmult/wiki/schubmult-home) for more details and examples.