# schubmult.rings.schubert_ring — Schubert ring classes and elements

Source: src/schubmult/rings/schubert_ring.py

Purpose
- Concrete ring classes exposing Schubert-basis algebras (single and double flavors) with dict-backed elements supporting arithmetic, pretty printing, and integrations with multiplication kernels.

Primary types
- DoubleSchubertRing — ring container that constructs DoubleSchubertElement instances and wires multiplication to kernels in mult.double.
- DoubleSchubertElement — dict-backed element representing Schubert-basis expressions; supports addition, multiplication, substitution, and canonical rendering.

Helpers and factories
- DSx — factory for basis elements.
- Integration points to elem_sym and polynomial basis utilities for coefficient expansions.

Notes
- Elements inherit from BaseSchubertElement; equality/approximate comparisons account for symbolic coefficient expansions and substitution rules.
