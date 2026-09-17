<a id="schubmult._scripts.diag_subwords"></a>

# schubmult.\_scripts.diag\_subwords

Compare per-subword weights between alphabet-vine (canonical_forest filter)
and reflection model (reduced+omega filter) on the same long_word index sets.

Both models select subwords of long_word(n). For each kept subword, print
(idx, values, per-letter weights, total weight) under each model.

