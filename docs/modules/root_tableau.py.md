# root_tableau.py



## _length_of_row(grid, row)



## _count_boxes(grid)



## _root_compare(root1, root2)



## _word_from_grid(grid0, as_grid, as_ordering, with_compatible_seq)

Two modes:
  - as_grid=True: return an object-array the same shape as grid0 where each
    occupied cell contains the recording letter (cell[1]) and empty cells
    are None.
  - as_grid=False: reconstruct the reduced word (sequence of simple-reflection
    indices) by repeatedly:
      * finding the occupied cell whose recording letter (cell[1]) is maximal
        and, among those, is farthest to the right (largest column index;
        break ties by largest row index),
      * popping that box, appending cell[0][0] to the collected word,
      * applying the root-shift corresponding to cell[0] to the region
        above the popped box and to the part of the same row left of the box,
      * repeating until no boxes remain.
    Returns a tuple(reversed(collected_letters)) to match the original reduced
    word orientation used elsewhere.

## _root_shift(root, spots)

Return a callable shift(grid_slice) -> new_grid_slice that applies the
appropriate root-reflection action to every non-None cell of the input
object-array slice. Uses Permutation.ref_product(...) as the reflection
provider and is defensive about act_root signatures.

## class RootTableau(CrystalGraph, GridPrint)

Root tableau with dual knuth equivalence

