<a id="schubmult._scripts.glide_lr_rule_nope"></a>

# schubmult.\_scripts.glide\_lr\_rule\_nope

<a id="schubmult._scripts.glide_lr_rule_nope.glide_rc_try"></a>

#### glide\_rc\_try

```python
def glide_rc_try(comp, beta, length)
```

GlidePoly polynomial of ``comp`` built by inverting the omega insertion.

We enumerate the compatible set-valued labelings ``kappa`` of the indexed
forest ``F = weak_composition_to_indfor(comp)`` (the glide definition), and
realize each labeling as a WCGraph by writing down its (word, compatible
sequence) pair *explicitly* -- the published inverse of the set-valued omega
insertion -- and feeding it to ``WCGraph.from_word_compatible``.

Concretely:

* The canonical left-binary-search labeling ``P`` of ``F`` and the decreasing
  labeling ``Q`` of the principal reduced RC graph form the omega pair for
  ``F``. The map ``Gamma = omega_reduced_word_from_labelings`` reads the
  reduced word ``W`` off ``(P, Q)`` -- this is the inverse of the insertion,
  so ``W`` is determined by ``F`` (not chosen at random). Node ``v`` sits at
  word position ``len(W) - Q(v)`` and carries the reduced letter
  ``ell(v) = W[len(W) - Q(v)]``.
* A labeling ``kappa`` places the letter ``ell(v)`` into every row
  ``r in kappa(v)``. Reading the resulting ``(row, letter)`` pairs in
  ``(row, -letter)`` order gives a weakly-increasing compatible sequence and
  a word whose rows are strictly decreasing, i.e. exactly the data
  ``WCGraph.from_word_compatible`` consumes. The multiset of compatible
  values is the disjoint union of the ``kappa(v)``, so the WCGraph monomial
  equals ``x^kappa``.

Each WCGraph contributes ``beta**(|kappa| - |F|)`` times its monomial; beta is
the degree ``-1`` homogenizer recording extra labels beyond one per node.

