<!-- filepath: src/schubmult/schub_lib/root_tableau.py -->

- `FunctionDef` — `_is_valid_outer_corner` (line 28)
    - Outer-corner predicate used by up_jdt_slide.
- `FunctionDef` — `_is_valid_inner_corner` (line 48)
    - Inner-corner predicate used by down_jdt_slide.
- `FunctionDef` — `_length_of_row` (line 66)
- `FunctionDef` — `_count_boxes` (line 70)
- `FunctionDef` — `_root_compare` (line 80)
- `FunctionDef` — `_word_from_grid` (line 92)
    - Two modes:
- `FunctionDef` — `_root_shift` (line 262)
    - Return a callable shift(grid_slice) -> new_grid_slice that applies the
- `FunctionDef` — `_validate_grid` (line 330)
    - Lightweight validation of a root-grid; raises on malformed cells.
- `FunctionDef` — `_snap_grid` (line 352)
    - Return a compact, JSON-like snapshot of the grid for debug messages.
- `ClassDef` — `RootTableau` (line 363)
    - Root tableau with dual knuth equivalence