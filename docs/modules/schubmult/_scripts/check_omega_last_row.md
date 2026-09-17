<a id="schubmult._scripts.check_omega_last_row"></a>

# schubmult.\_scripts.check\_omega\_last\_row

<a id="schubmult._scripts.check_omega_last_row.check_disjoint_vs_squash_omega"></a>

#### check\_disjoint\_vs\_squash\_omega

```python
def check_disjoint_vs_squash_omega(n: int) -> tuple[bool, int, int]
```

Compare omega_invariant[1] for disjoint_union vs squash_product.

Loop over minimal-length RC graphs rc. For each rc, loop over elementary
symmetric RC graphs of matching length and last descent, and verify:

    rc.disjoint_union(elem_sym_rc).omega_invariant[1]
    ==
    rc.squash_product(elem_sym_rc).omega_invariant[1]

