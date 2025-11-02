# schubmult.schub_lib.perm_lib — permutation primitives

Source: src/schubmult/schub_lib/perm_lib.py

Overview
- Lightweight Permutation class and utilities for Lehmer/code conversions, Bruhat relations, cycles, and coset representatives used across the codebase.

Primary symbol
- Permutation — core permutation type with helpers:
  - construction: from code/cycles, uncode, code
  - properties: inv, code/trimcode, descents, cycles
  - comparators and Bruhat utilities used by algorithms that enumerate or compare permutations.

Useful utilities
- code, uncode, permtrim, theta, strict_theta, longest_element — common transformations.
- cyclic_sort, cyclic_sort_min, get_cycles — cycle and ordering helpers.

Notes
- Permutation objects are treated as immutable and often cached for allocation savings; many algorithms rely on fast equality and code-based comparisons.
