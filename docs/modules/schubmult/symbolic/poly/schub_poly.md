<a id="schubmult.symbolic.poly.schub_poly"></a>

# schubmult.symbolic.poly.schub\_poly

<a id="schubmult.symbolic.poly.schub_poly.isobaric_strip_on_dschub_dict"></a>

#### isobaric\_strip\_on\_dschub\_dict

```python
def isobaric_strip_on_dschub_dict(start, length, perm_dict, coeff_genset,
                                  beta)
```

Apply one isobaric strip to a whole ``{perm: coeff}`` dict, folded.

Coefficients landing on the same permutation merge at every stage instead of
being carried per input basis element, mirroring ``compute_vpathdicts``.

<a id="schubmult.symbolic.poly.schub_poly.apply_isobaric_to_schub_dict"></a>

#### apply\_isobaric\_to\_schub\_dict

```python
def apply_isobaric_to_schub_dict(diff_perm, perm_dict, coeff_genset, beta)
```

Fold every strip of ``diff_perm`` over the whole dict, merging between strips.

<a id="schubmult.symbolic.poly.schub_poly.groth_elem_as_schub_dict"></a>

#### groth\_elem\_as\_schub\_dict

```python
@cache
def groth_elem_as_schub_dict(perm, beta)
```

Expand a Grothendieck basis element G_perm into Schubert basis terms.

Uses the strip-isobaric construction from the groth_elem_as_schub method.
Returns ``{Permutation: coeff}``.

