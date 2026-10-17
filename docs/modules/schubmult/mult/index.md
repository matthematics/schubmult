<a id="schubmult.mult"></a>

# schubmult.mult

Multiplication algorithms for Schubert polynomials.

This module provides kernels for computing products of Schubert polynomials
in various settings:

- schubmult_py: Ordinary (single) Schubert polynomial multiplication
- schubmult_double: Double Schubert polynomial multiplication
- schubmult_q: Quantum Schubert polynomial multiplication
- schubmult_q_double: Quantum double Schubert polynomial multiplication
- grothmult_double: Double Grothendieck multiplication by a degree-one class
- grothmult_q_double: Quantum double Grothendieck multiplication (conjectural Molev--Sagan rule)
- grothmult_q: Quantum (single) Grothendieck multiplication

Also includes positivity utilities (posify, compute_positive_rep) for root-based representations.
Exports resolve lazily (PEP 562) so that importing one kernel does not load all the others.

