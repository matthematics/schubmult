# scripts.schubmult_double — CLI for double Schubert products

Source: src/scripts/schubmult_double.py

Purpose
- Command-line tool for computing and pretty-printing double Schubert polynomial products; supports display-positive mode via posify.

Important functions
- _display_full — format and print full coefficient dictionaries; prepare substitution maps for nicer variable rendering.
- pre_posify / sv_posify — wrappers that convert raw symbolic coefficients into display-positive canonical forms by delegating to mult.positivity.posify.
- main — CLI entrypoint wired by pyproject scripts (schubmult_double).

Usage notes
- The script supports flags for coproducts, display-positive output, and variable control. See tests in src/tests/script_tests for example invocations.
