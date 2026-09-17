<a id="schubmult._scripts.count_sanity"></a>

# schubmult.\_scripts.count\_sanity

Count subwords of long_word(n) under two filters and compare.

<a id="schubmult._scripts.count_sanity.count_sylvester_subwords"></a>

#### count\_sylvester\_subwords

```python
def count_sylvester_subwords(code, n=None)
```

Count `1`: subwords of full (barred+unbarred) long word in Syl(F).

<a id="schubmult._scripts.count_sanity.count_reduced_with_forest"></a>

#### count\_reduced\_with\_forest

```python
def count_reduced_with_forest(code, n=None)
```

Count `2`: subwords of the FULL long word (barred+unbarred) whose value
sequence is a reduced word for perm = uncode(code) and whose
omega_invariant[0] (of reversed values) matches the principal-RC target.

