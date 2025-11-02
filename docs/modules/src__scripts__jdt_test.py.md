# src/scripts/jdt_test.py

Test whether raising_operator / lowering_operator commute with
(arbitrary) sequences of up_jdt_slide followed by rectification.

Usage: run this script from the repo root (with the correct PYTHONPATH / venv
active) to run randomized trials and report any counterexamples.

## apply_up_seq_and_rect(rt, seq)

For a given RootTableau (assumed straight), apply a sequence of "create hole
at (i,j) then up_jdt_slide(i,j)" operations and finally rectify().

Each (i,j) may extend the tableau (a hole on an outer corner). We create
that hole by extending the underlying grid and setting that cell to None,
then calling up_jdt_slide on the new tableau.

## grids_equal(a, b)

Compare two RootTableau by their underlying root grids (shape + cell equality).

## test_one_case(T, index, op_name, rc)

Test commutation for one operator index:
 left = op( Rect( UpSeq(T) ) )
 right = Rect( UpSeq( op(T) ) )
Returns (passed, message)

Note: do not suppress exceptions from raising_operator / lowering_operator —
these should return None when the operator is not defined.

## random_up_seq(rt, max_len)

Generate a valid sequence of up-jdt hole positions for tableau `rt`.

Each chosen (i,j) satisfies the outer-corner predicate above w.r.t. the
current tableau. After choosing a hole we perform the up_jdt_slide to
update the tableau so subsequent choices remain valid.

## run_random_tests(num_cases)

Run randomized tests and log the outcome of each trial (both raise and lower).
Ensure the up-jdt sequence is actually valid for the sampled tableau: regenerate
until apply_up_seq_and_rect(T, seq) succeeds (or give up after attempts).

## run_complete_tests()

Run randomized tests and log the outcome of each trial (both raise and lower).
Ensure the up-jdt sequence is actually valid for the sampled tableau: regenerate
until apply_up_seq_and_rect(T, seq) succeeds (or give up after attempts).

