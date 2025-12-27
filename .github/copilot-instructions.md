# AI Coding Instructions for schubmult

## Project Purpose
`schubmult` computes Littlewood-Richardson coefficients for Schubert polynomials (ordinary, double, quantum, quantum double) with optional parabolic subgroup support. The package integrates with SymPy and provides CLI tools for algebraic computation in Schubert calculus.

## Architecture Overview

### Core Components
1. **Combinatorial Structures** ([src/schubmult/schub_lib/](src/schubmult/schub_lib/))
   - `Permutation`: Central object indexed by permutations, uses Lehmer codes (`code`/`trimcode` attributes)
   - `RCGraph`: Reduced-compatible graphs with crystal structure (raising/lowering operators)
   - `RootTableau`: Dual Knuth equivalence implementation with JDT slides and Edelman-Greene invariants
   - `CrystalGraph`: Abstract base for crystal operators (`raising_operator(i)`, `lowering_operator(i)`)

2. **Ring Structures** ([src/schubmult/rings/](src/schubmult/rings/))
   - `BaseSchubertRing`: Parent ring for Schubert element types
   - `DoubleSchubertRing`, `SingleSchubertRing`: Concrete implementations for single/double Schubert polynomials
   - `QuantumDoubleSchubertRing`: Quantum variants with q-variable bookkeeping
   - Elements stored as dicts `{Permutation: coefficient}` with ring operations

3. **Multiplication Kernels** ([src/schubmult/mult/](src/schubmult/mult/))
   - `schubmult_py`: Ordinary Schubert polynomial multiplication via theta codes and v-path dictionaries
   - `schubmult_double`: Double Schubert polynomials with y/z variable support
   - `schubmult_q`: Quantum multiplication with parabolic support (Peterson-Woodward theorem)
   - `schubmult_q_double`: Quantum double (conjectural for most cases)

4. **CLI Scripts** ([src/schubmult/_scripts/](src/schubmult/_scripts/))
   - All use `schub_argparse()` for consistent argument parsing
   - Entry points: `schubmult_py`, `schubmult_double`, `schubmult_q`, `schubmult_q_double`, `lr_rule_verify`
   - Permutations specified as space-separated integers or via `--code` flag for Lehmer codes

### Critical Data Flows
- **Permutation → RCGraph → RootTableau**: Main conversion chain for combinatorial operations
- **Crystal operators**: Apply raising/lowering ops on crystal structures, preserve EG invariants
- **Ring multiplication**: Dict-based accumulation using `add_perm_dict()` utility
- **Positivity display**: Integer linear programming via PuLP for `--display-positive` mode

## Development Workflows

### Testing
```bash
pytest                           # Run all tests
pytest tests/ring_tests/         # Ring-specific tests
pytest tests/schub_tests/        # Combinatorial structure tests
pytest -k "test_root_tableau"    # Specific test pattern
```
- Test config in `pytest.ini` and `[tool.pytest.ini_options]` (pyproject.toml)
- Uses `--import-mode=importlib` to handle src-layout
- Script tests in `tests/script_tests/` with JSON test case files

### Building & Installation
```bash
pip install -e .                 # Editable dev install
python -m build                  # Build distribution packages
```
- Uses `setuptools>=61.0` backend with PEP 621 metadata in `pyproject.toml`
- Dynamic versioning from `schubmult.__version__` (currently `"4.0.0dev"`)
- CLI scripts auto-registered via `[project.scripts]` section

### Running Scripts
```bash
schubmult_py 3 1 2 - 2 1 3                    # Ordinary Schubert product
schubmult_double 3 1 2 - 2 1 3 --display-positive  # Double with positivity
schubmult_q 3 1 2 - 2 1 3 --parabolic 1       # Quantum with parabolic
schubmult_py --code 2 0 - 1 0                 # Using Lehmer codes
```
- All scripts support `--code` for Lehmer code input, `--no-print` to suppress output
- Double/quantum variants have `--display-positive` for root-based display (uses MILP)
- Quantum scripts accept `--parabolic` for parabolic subgroup generators

## Project-Specific Conventions

### Symbolic Computation Layer
- **Hybrid SymPy/SymEngine**: Uses `symengine` for fast arithmetic, `sympy` for display/latex
- Import pattern: `from schubmult.symbolic import sympify, Add, Mul, S` (wraps both libraries)
- Variable generation: `GeneratingSet("x")` creates indexed variable sets (x₁, x₂, ...)
- DO NOT import directly from sympy/symengine; use `schubmult.symbolic` facade

### Permutation Representations
- **Array form**: `[3, 1, 2]` means π(1)=3, π(2)=1, π(3)=2
- **Lehmer code**: `perm.code` returns full code, `perm.trimcode` strips trailing zeros
- Conversion: `uncode(lehmer_list)` → `Permutation`, `Permutation([3,1,2]).code` → code
- Descents: `perm.descents()` returns list of descent positions

### Ring Element Operations
- Elements are dict subclasses: `{Permutation([2,1,3]): coeff, ...}`
- Multiplication via `*` operator dispatches to multiplication kernels
- Expansion: `.expand()` converts to explicit polynomial in variables
- Utility: `add_perm_dict(dict1, dict2)` merges coefficient dictionaries

### Crystal Structure Patterns
- Crystal operators return `None` if undefined (not an error condition)
- Always check `if result is None:` before using operator results
- Highest/lowest weight: `.to_highest_weight()` returns `(hw_element, sequence)`
- Weight vectors: `.crystal_weight` property (tuple of integers)

### Code Style
- Line length: 200 chars (configured in ruff)
- Imports: Use absolute imports `from schubmult.X import Y`
- Lazy imports in `__init__.py`: Uses `__getattr__` for heavy dependencies
- Ignore test files for linting (see `[tool.ruff.lint.per-file-ignores]`)

## Common Pitfalls

1. **Crystal operator validation**: Always verify `.edelman_greene_invariant` is preserved in root tableau operations
2. **Permutation indexing**: 1-indexed for descent positions but 0-indexed in Python arrays
3. **RC-graph normalization**: Use `.normalize()` after manual construction
4. **Quantum q-variable**: Stored in separate coefficient ring, factor via `factor_out_q()`
5. **MILP display**: `--display-positive` requires PuLP and can be slow; use `--optimizer-message` for debug

## Key File References

- Argument parsing: [src/schubmult/utils/argparse.py](src/schubmult/utils/argparse.py) (`schub_argparse()`)
- Permutation utilities: [src/schubmult/schub_lib/perm_lib.py](src/schubmult/schub_lib/perm_lib.py)
- RC-graph operations: [src/schubmult/schub_lib/rc_graph.py](src/schubmult/schub_lib/rc_graph.py) (1419 lines, see `monk_crystal_mul`)
- Root tableau logic: [src/schubmult/schub_lib/root_tableau.py](src/schubmult/schub_lib/root_tableau.py) (EG invariant tracking)
- Main multiplication: [src/schubmult/mult/single.py](src/schubmult/mult/single.py) (`schubmult_py` kernel)

## External Dependencies

- **Core**: numpy, sympy>=1.14, symengine>=0.14, cachetools, sortedcontainers
- **Optimization**: PuLP>=2.7 (for positivity display via MILP)
- **Optional**: sagemath-standard (for Sage integration, `[sage]` extra)
- **Build**: setuptools>=61, joblib (parallelization), psutil (resource monitoring)

## Verification Workflows

The `lr_rule_verify` script runs exhaustive LR coefficient verification:
- Generates all permutations up to size n
- Compares algorithmic results against classical formulas
- Outputs JSON with verification records (see `logs/*.verification.json`)
- Uses multiprocessing (specify `num_processors` argument)
- Resumable: Loads existing `.verification.json` to continue interrupted runs
