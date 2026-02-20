# schubmult

## Program and package for rapid computation of Littlewood-Richardson coefficients of Schubert polynomials, compliant with sympy/symengine (and hence indirectly Sage)

The main purpose of this python package is for doing Schubert calculus-related calculations in python (and/or Sage). 

- Kinds of things covered (not exhaustive):
    - Permutation library
	- Fast multiplication of single, double, mixed-variable Schubert polynomials; quantum, quantum double, quantum mixed-variable, and all parabolic versions.
	- Noncommutative algebras such as NSym and the free algebra on words of nonnegative integers augmented with combinatorial bases.
	- RC graphs/PDs, BPDs, HPDs, SSYT, EG tableaux, and algebraic structures derived from them (Coxeter-Knuth insertion, RSK, RC graph transition formulas, tableaux decompositions). Kashiwara/Demazure crystal raising/lowering operators
    - Compatible with sympy and symengine, and hence Sage, probably not terribly difficult to integrate with libraries I'm not familiar with.


[Docs to be hosted on Wiki](https://github.com/matthematics/schubmult/wiki/schubmult-home)


## To install dev version

```
pip install git+https://github.com/matthematics/schubmult.git
```

## RCGraph and BPD Functionality

The package implements two main combinatorial models for Schubert calculus:

- **RCGraph (Reduced Compatible Graphs):** Encodes reduced words for permutations as graphs, supporting crystal operations and algebraic manipulations.
- **BPD (Bumpless Pipe Dreams):** Represents tilings of an $n \times n$ grid with local tile rules, providing an alternative model for Schubert polynomials.
- **HPD (Hybrid Pipe Dreams):** They exist in the library and are functional, not nearly as well developed at this time.

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
- RCGraphs have a product through RCGraphRing (similar to concatenation, not polynomial product)

### Example Usage

```python
from schubmult import RCGraph, BPD, Permutation

rc = RCGraph.random_rc_graph(Permutation([5,1,6,2,4,3]), 5)
bpd = BPD.from_rc_graph(rc)
rc2 = bpd.to_rc_graph()
assert rc2 == rc
```
