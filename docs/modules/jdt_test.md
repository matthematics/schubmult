<!-- filepath: src/scripts/jdt_test.py -->

Test whether raising_operator / lowering_operator commute with
(arbitrary) sequences of up_jdt_slide followed by rectification.

Usage: run this script from the repo root (with the correct PYTHONPATH / venv
active) to run randomized trials and report any counterexamples.

- `FunctionDef` — `apply_up_seq_and_rect` (line 35)
    - For a given RootTableau (assumed straight), apply a sequence of "create hole
- `FunctionDef` — `grids_equal` (line 55)
    - Compare two RootTableau by their underlying root grids (shape + cell equality).
- `FunctionDef` — `test_one_case` (line 94)
    - Test commutation for one operator index:
- `FunctionDef` — `random_up_seq` (line 137)
    - Generate a valid sequence of up-jdt hole positions for tableau `rt`.
- `FunctionDef` — `run_random_tests` (line 161)
    - Run randomized tests and log the outcome of each trial (both raise and lower).
- `FunctionDef` — `run_complete_tests` (line 207)
    - Run randomized tests and log the outcome of each trial (both raise and lower).