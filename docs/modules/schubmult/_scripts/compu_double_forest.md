<a id="schubmult._scripts.compu_double_forest"></a>

# schubmult.\_scripts.compu\_double\_forest

<a id="schubmult._scripts.compu_double_forest.class_key"></a>

#### class\_key

```python
def class_key(values)
```

Omega-class invariant: P-symbol of omega_insertion.

<a id="schubmult._scripts.compu_double_forest.build_match"></a>

#### build\_match

```python
def build_match(comp, n)
```

Match B-subwords to A-subwords by omega_insertion Q-symbol.

A := Sylvester subwords of long_word (canonical_forest == sylv_canon).
B := reduced-word subwords for uncode(comp) in the omega-class of the
     canonical reduced word of perm.

For each subword W (A or B), compute (P,Q) = omega_insertion(reverse(W)).
Both A and B subwords have shape F. Two subwords with identical Q place
their i-th reversed letter at the same node, so a_rev[i] <-> b_rev[i]
pairs them position-by-position; equivalently a[pos] <-> b[pos] in
original (long_word) order. The Q-bijection is unique and total.

