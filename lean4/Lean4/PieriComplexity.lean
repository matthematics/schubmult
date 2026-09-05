/-
Copyright (c) 2026 Matthew J. Samuel. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Matthew J. Samuel
-/
import Mathlib

/-!
# Complexity of the `pieri` algorithm

Formalization of Proposition `prop:piericomplexity` ("Complexity of **pieri**") from
`writing/schubmult.tex`.  The proposition reads: with `M := |T|` the number of admissible
`k`-chains from `u` of length at most `p` returned by `pieri u p k`, the algorithm runs in
`O(k (n-k) n · M)` elementary operations, that is `O(n³)` per enumerated permutation, and
`M ≤ #{w ∈ Sₙ | w ≥ u, ℓ(w) ≤ ℓ(u) + p}`.

Both halves are proved, `pieri_complexity` being the combined statement.  The only thing
assumed is multiplicity-freeness (Lemma `lemma:pierimultfree`), quoted from the literature
and recorded as the axiom `table_nodup`; `#print axioms` reports it.

## The cost model

`stepCost` charges one unit per element visited by each of the three nested traversals of
`stepNode`, which is written to follow Step 2 of Method `method:pieri` literally: scan the
`k` left positions to build `A`, scan the `n - k` right positions for `b`, and for each
surviving pair scan the positions between `a` and `b` for the Bruhat-cover test.  Each
summand of `stepCost` is therefore indexed by a list that `stepNode` actually traverses,
and `stepCost_le` recovers the paper's `k + (n-k) · k · n` shape from those list lengths.

## Length and Bruhat order

Neither is in Mathlib, so both are built here.  `len` is the inversion count, the Coxeter
length of `Sₙ`.  `BruhatLE` is the transposition characterization of the Bruhat order:
`u ≤ w` when `w` is reached from `u` by transpositions each raising `len` by one.  The
combinatorial content is `len_rswap`: a `k`-strip cover raises the inversion count by
exactly one.  This is what puts the endpoints in the interval the proposition names.

## Indexing

Positions and values are `0`-indexed, so the paper's `1 ≤ a ≤ k < b ≤ n` reads `a < k ≤ b`,
and the paper's `β = ∞` is the sentinel bound `n + 1`, which exceeds every value `w b < n`.
-/

namespace Schubmult.Pieri

variable {n : ℕ}

/-- Permutations of `{1, …, n}`, the ambient group `Sₙ`. -/
abbrev Perm (n : ℕ) := Equiv.Perm (Fin n)

/-- `w t_{ab}`: the transposition of *positions* `a` and `b`, acting on the right. -/
def rswap (w : Perm n) (a b : Fin n) : Perm n := w * Equiv.swap a b

theorem rswap_apply (w : Perm n) (a b x : Fin n) :
    rswap w a b x = w (Equiv.swap a b x) := rfl

theorem rswap_left (w : Perm n) (a b : Fin n) : rswap w a b a = w b := by
  simp [rswap_apply]

theorem rswap_right (w : Perm n) (a b : Fin n) : rswap w a b b = w a := by
  simp [rswap_apply]

theorem rswap_other (w : Perm n) {a b x : Fin n} (hx : x ≠ a) (hx' : x ≠ b) :
    rswap w a b x = w x := by
  simp [rswap_apply, Equiv.swap_apply_of_ne_of_ne hx hx']

/-! ### Length

`len` is the number of inversions, the Coxeter length of `Sₙ`. -/

/-- The inversion set of `w`. -/
def invs (w : Perm n) : Finset (Fin n × Fin n) :=
  Finset.univ.filter fun q => q.1 < q.2 ∧ w q.2 < w q.1

/-- `ℓ(w)`, the number of inversions of `w`. -/
def len (w : Perm n) : ℕ := (invs w).card

theorem mem_invs {w : Perm n} {q : Fin n × Fin n} :
    q ∈ invs w ↔ q.1 < q.2 ∧ w q.2 < w q.1 := by
  simp [invs]

/-! ### Bruhat covers -/

/-- The positions strictly between `a` and `b`, in the order the cover test scans them. -/
def between (a b : Fin n) : List (Fin n) :=
  (List.finRange n).filter fun c => a < c ∧ c < b

theorem mem_between {a b c : Fin n} : c ∈ between a b ↔ a < c ∧ c < b := by
  simp [between]

theorem length_between_le (a b : Fin n) : (between a b).length ≤ n := by
  have h := List.length_filter_le (fun c : Fin n => decide (a < c ∧ c < b)) (List.finRange n)
  simpa [between, List.length_finRange] using h

/-- `w ⋖ w t_{ab}` in Bruhat order: the values at `a` and `b` are in order, and no position
strictly between them carries a value strictly between them. -/
def IsCover (w : Perm n) (a b : Fin n) : Prop :=
  w a < w b ∧ ∀ c ∈ between a b, ¬(w a < w c ∧ w c < w b)

instance (w : Perm n) (a b : Fin n) : Decidable (IsCover w a b) := by
  unfold IsCover; infer_instance

/-! ### A cover raises the inversion count by exactly one

Write `σ = t_{ab}`.  On pairs `(i, j)` with `i < j` that `σ` does not reverse, transporting
by `σ` matches inversions of `w σ` with inversions of `w`.  The pairs that `σ` does reverse
are exactly `(a, b)` and, for each `c` with `a < c < b`, the two pairs `(a, c)` and `(c, b)`;
on those the cover hypothesis says `w` and `w σ` have *the same* inversions, so `(a, b)` is
the single new one. -/

/-- Transport of a pair by `σ = t_{ab}`, sorted; the identity on the pairs `σ` reverses. -/
private def tw (a b : Fin n) (q : Fin n × Fin n) : Fin n × Fin n :=
  if Equiv.swap a b q.2 < Equiv.swap a b q.1 then q
  else (Equiv.swap a b q.1, Equiv.swap a b q.2)

/-- The pairs `i < j` that `t_{ab}` reverses are `(a, b)` and the two families straddling a
position strictly between `a` and `b`. -/
private theorem rev_cases {a b i j : Fin n} (hab : a < b) (hij : i < j)
    (h : Equiv.swap a b j < Equiv.swap a b i) :
    (i = a ∧ j = b) ∨ (i = a ∧ a < j ∧ j < b) ∨ (j = b ∧ a < i ∧ i < b) := by
  rcases eq_or_ne i a with rfl | hi
  · rcases eq_or_ne j b with rfl | hj
    · exact Or.inl ⟨rfl, rfl⟩
    · rw [Equiv.swap_apply_left, Equiv.swap_apply_of_ne_of_ne hij.ne' hj] at h
      exact Or.inr (Or.inl ⟨rfl, hij, h⟩)
  · rcases eq_or_ne i b with rfl | hib
    · rw [Equiv.swap_apply_right,
        Equiv.swap_apply_of_ne_of_ne (hab.trans hij).ne' hij.ne'] at h
      exact absurd (hab.trans hij) (asymm h)
    · rw [Equiv.swap_apply_of_ne_of_ne hi hib] at h
      rcases eq_or_ne j a with rfl | hja
      · rw [Equiv.swap_apply_left] at h
        exact absurd (hij.trans hab) (asymm h)
      · rcases eq_or_ne j b with rfl | hjb
        · rw [Equiv.swap_apply_right] at h
          exact Or.inr (Or.inr ⟨rfl, h, hij⟩)
        · rw [Equiv.swap_apply_of_ne_of_ne hja hjb] at h
          exact absurd hij (asymm h)

/-- On a reversed pair, `w t_{ab}` has the same inversions as `w`, except that `(a, b)`
becomes one.  This is where the cover hypothesis is used. -/
private theorem mem_invs_rswap_rev (w : Perm n) {a b : Fin n} (hab : a < b)
    (hc : IsCover w a b) {q : Fin n × Fin n} (h1 : q.1 < q.2)
    (h2 : Equiv.swap a b q.2 < Equiv.swap a b q.1) :
    q ∈ invs (rswap w a b) ↔ (q ∈ invs w ∨ q = (a, b)) := by
  obtain ⟨i, j⟩ := q
  simp only at h1 h2
  rcases rev_cases hab h1 h2 with ⟨rfl, rfl⟩ | ⟨rfl, haj, hjb⟩ | ⟨rfl, hai, hib⟩
  · simp [mem_invs, rswap_left, rswap_right, hab, hc.1]
  · have hno := hc.2 j (mem_between.mpr ⟨haj, hjb⟩)
    have hjne : j ≠ i := haj.ne'
    have hjne' : j ≠ b := hjb.ne
    simp only [mem_invs, rswap_left, rswap_other w hjne hjne', Prod.mk.injEq]
    constructor
    · rintro ⟨-, hlt⟩
      refine Or.inl ⟨h1, ?_⟩
      rcases lt_trichotomy (w j) (w i) with h | h | h
      · exact h
      · exact absurd (w.injective h) hjne
      · exact absurd ⟨h, hlt⟩ hno
    · rintro (⟨-, hlt⟩ | ⟨-, hjb'⟩)
      · exact ⟨h1, hlt.trans hc.1⟩
      · exact absurd hjb' hjne'
  · have hno := hc.2 i (mem_between.mpr ⟨hai, hib⟩)
    have hine : i ≠ a := hai.ne'
    have hine' : i ≠ j := hib.ne
    simp only [mem_invs, rswap_right, rswap_other w hine hine', Prod.mk.injEq]
    constructor
    · rintro ⟨-, hlt⟩
      refine Or.inl ⟨h1, ?_⟩
      rcases lt_trichotomy (w j) (w i) with h | h | h
      · exact h
      · exact absurd (w.injective h) (Ne.symm hine')
      · exact absurd ⟨hlt, h⟩ hno
    · rintro (⟨-, hlt⟩ | ⟨hia, -⟩)
      · exact ⟨h1, hc.1.trans hlt⟩
      · exact absurd hia hine

/-- On a pair `t_{ab}` does not reverse, transport by `t_{ab}` matches the inversion
sets. -/
private theorem mem_invs_rswap_nonrev (w : Perm n) {a b : Fin n} {q : Fin n × Fin n}
    (h1 : q.1 < q.2) (h2 : ¬(Equiv.swap a b q.2 < Equiv.swap a b q.1)) :
    q ∈ invs (rswap w a b) ↔ (Equiv.swap a b q.1, Equiv.swap a b q.2) ∈ invs w := by
  have hne : Equiv.swap a b q.1 ≠ Equiv.swap a b q.2 := fun h =>
    absurd ((Equiv.swap a b).injective h) h1.ne
  have hlt : Equiv.swap a b q.1 < Equiv.swap a b q.2 := lt_of_le_of_ne (not_lt.mp h2) hne
  simp only [mem_invs, rswap_apply, hlt, h1, true_and]

private theorem tw_rev {a b : Fin n} {q : Fin n × Fin n}
    (h : Equiv.swap a b q.2 < Equiv.swap a b q.1) : tw a b q = q := by
  simp [tw, h]

private theorem tw_nonrev {a b : Fin n} {q : Fin n × Fin n}
    (h : ¬(Equiv.swap a b q.2 < Equiv.swap a b q.1)) :
    tw a b q = (Equiv.swap a b q.1, Equiv.swap a b q.2) := by
  simp [tw, h]

theorem invs_rswap (w : Perm n) {a b : Fin n} (hab : a < b) (hc : IsCover w a b) :
    invs (rswap w a b) = insert (a, b) ((invs w).image (tw a b)) := by
  classical
  ext q
  simp only [Finset.mem_insert, Finset.mem_image]
  constructor
  · intro hq
    have h1 : q.1 < q.2 := (mem_invs.mp hq).1
    by_cases h2 : Equiv.swap a b q.2 < Equiv.swap a b q.1
    · rcases (mem_invs_rswap_rev w hab hc h1 h2).mp hq with hin | rfl
      · exact Or.inr ⟨q, hin, tw_rev h2⟩
      · exact Or.inl rfl
    · refine Or.inr ⟨(Equiv.swap a b q.1, Equiv.swap a b q.2),
        (mem_invs_rswap_nonrev w h1 h2).mp hq, ?_⟩
      rw [tw_nonrev (by simpa [Equiv.swap_apply_self] using asymm h1)]
      simp [Equiv.swap_apply_self]
  · rintro (rfl | ⟨r, hr, rfl⟩)
    · exact mem_invs.mpr ⟨hab, by rw [rswap_left, rswap_right]; exact hc.1⟩
    · have h1 : r.1 < r.2 := (mem_invs.mp hr).1
      by_cases h2 : Equiv.swap a b r.2 < Equiv.swap a b r.1
      · rw [tw_rev h2]
        exact (mem_invs_rswap_rev w hab hc h1 h2).mpr (Or.inl hr)
      · rw [tw_nonrev h2]
        refine (mem_invs_rswap_nonrev w (a := a) (b := b) ?_ ?_).mpr ?_
        · exact lt_of_le_of_ne (not_lt.mp h2)
            fun h => absurd ((Equiv.swap a b).injective h) h1.ne
        · simpa [Equiv.swap_apply_self] using asymm h1
        · simpa [Equiv.swap_apply_self] using hr

theorem len_rswap (w : Perm n) {a b : Fin n} (hab : a < b) (hc : IsCover w a b) :
    len (rswap w a b) = len w + 1 := by
  classical
  have hinj : Set.InjOn (tw a b) (invs w) := by
    intro r hr s hs h
    have hr1 : r.1 < r.2 := (mem_invs.mp hr).1
    have hs1 : s.1 < s.2 := (mem_invs.mp hs).1
    by_cases hR : Equiv.swap a b r.2 < Equiv.swap a b r.1 <;>
      by_cases hS : Equiv.swap a b s.2 < Equiv.swap a b s.1
    · rwa [tw_rev hR, tw_rev hS] at h
    · rw [tw_rev hR, tw_nonrev hS] at h
      rw [h] at hR
      simp only [Equiv.swap_apply_self] at hR
      exact absurd hR (asymm hs1)
    · rw [tw_rev hS, tw_nonrev hR] at h
      rw [← h] at hS
      simp only [Equiv.swap_apply_self] at hS
      exact absurd hS (asymm hr1)
    · rw [tw_nonrev hR, tw_nonrev hS] at h
      have h1 := congrArg Prod.fst h
      have h2 := congrArg Prod.snd h
      simp only at h1 h2
      exact Prod.ext ((Equiv.swap a b).injective h1) ((Equiv.swap a b).injective h2)
  have hnot : (a, b) ∉ (invs w).image (tw a b) := by
    simp only [Finset.mem_image, not_exists, not_and]
    intro r hr hEq
    have hr1 : r.1 < r.2 := (mem_invs.mp hr).1
    by_cases hR : Equiv.swap a b r.2 < Equiv.swap a b r.1
    · rw [tw_rev hR] at hEq
      subst hEq
      exact absurd (mem_invs.mp hr).2 (asymm hc.1)
    · rw [tw_nonrev hR] at hEq
      have h1 := congrArg Prod.fst hEq
      have h2 := congrArg Prod.snd hEq
      simp only at h1 h2
      have e1 : r.1 = b := by
        have := congrArg (Equiv.swap a b) h1
        simpa [Equiv.swap_apply_self] using this
      have e2 : r.2 = a := by
        have := congrArg (Equiv.swap a b) h2
        simpa [Equiv.swap_apply_self] using this
      rw [e1, e2] at hr1
      exact absurd hr1 (asymm hab)
  rw [len, invs_rswap w hab hc, Finset.card_insert_of_notMem hnot,
    Finset.card_image_of_injOn hinj, len]

/-- Bruhat order, in its transposition characterization: `u ≤ w` when `w` is obtained from
`u` by successively multiplying by transpositions, each step raising the length by one. -/
inductive BruhatLE {n : ℕ} : Perm n → Perm n → Prop
  | refl (w : Perm n) : BruhatLE w w
  | step {u w : Perm n} {a b : Fin n} :
      BruhatLE u w → len (rswap w a b) = len w + 1 → BruhatLE u (rswap w a b)

theorem bruhat_trans {u v w : Perm n} (h₁ : BruhatLE u v) (h₂ : BruhatLE v w) :
    BruhatLE u w := by
  induction h₂ with
  | refl => exact h₁
  | step _ hlen ih => exact BruhatLE.step ih hlen

/-! ### The algorithm

`stepNode` follows Step 2 of Method `method:pieri` literally: build `A`, loop over `b`, and
for each pair run the Bruhat-cover test. -/

/-- The `k` left positions scanned when building `A` (Step 2(a)). -/
def scanA (n k : ℕ) : List (Fin n) := (List.finRange n).filter fun a => (a : ℕ) < k

/-- The `n - k` right positions scanned for `b` (Step 2(b)). -/
def scanB (n k : ℕ) : List (Fin n) := (List.finRange n).filter fun b => k ≤ (b : ℕ)

theorem mem_scanA {k : ℕ} {a : Fin n} (h : a ∈ scanA n k) : (a : ℕ) < k := by
  simpa [scanA] using h

theorem mem_scanB {k : ℕ} {b : Fin n} (h : b ∈ scanB n k) : k ≤ (b : ℕ) := by
  simpa [scanB] using h

theorem length_scanA_le (k : ℕ) : (scanA n k).length ≤ k := by
  classical
  have hnd : (scanA n k).Nodup := (List.nodup_finRange n).filter _
  have h1 : (scanA n k).length = (scanA n k).toFinset.card :=
    (List.toFinset_card_of_nodup hnd).symm
  have h2 : (scanA n k).toFinset.card ≤ (Finset.range k).card := by
    refine Finset.card_le_card_of_injOn Fin.val (fun a ha => ?_)
      (fun x _ y _ h => Fin.val_injective h)
    exact Finset.mem_range.mpr (mem_scanA (List.mem_toFinset.mp ha))
  rw [h1]
  simpa using h2

theorem length_scanB_le (k : ℕ) : (scanB n k).length ≤ n - k := by
  classical
  have hnd : (scanB n k).Nodup := (List.nodup_finRange n).filter _
  have h1 : (scanB n k).length = (scanB n k).toFinset.card :=
    (List.toFinset_card_of_nodup hnd).symm
  have h2 : (scanB n k).toFinset.card ≤ (Finset.Ico k n).card := by
    refine Finset.card_le_card_of_injOn Fin.val (fun b hb => ?_)
      (fun x _ y _ h => Fin.val_injective h)
    exact Finset.mem_Ico.mpr ⟨mem_scanB (List.mem_toFinset.mp hb), b.isLt⟩
  rw [h1]
  simpa [Nat.card_Ico] using h2

/-- A frontier entry: an intermediate permutation together with its current bound. -/
structure Node (n : ℕ) where
  perm : Perm n
  bound : ℕ
  deriving DecidableEq

/-- `A`, the eligible left positions of Step 2(a). -/
def eligible (k : ℕ) (nd : Node n) : List (Fin n) :=
  (scanA n k).filter fun a => (nd.perm a : ℕ) < nd.bound

theorem length_eligible_le (k : ℕ) (nd : Node n) :
    (eligible k nd).length ≤ (scanA n k).length :=
  List.length_filter_le _ _

theorem mem_eligible {k : ℕ} {nd : Node n} {a : Fin n} (h : a ∈ eligible k nd) :
    a ∈ scanA n k :=
  List.mem_of_mem_filter h

/-- Step 2 of Method `method:pieri` applied to one frontier node. -/
def stepNode (k : ℕ) (nd : Node n) : List (Node n) :=
  (scanB n k).flatMap fun b =>
    if (nd.perm b : ℕ) < nd.bound then
      (eligible k nd).filterMap fun a =>
        if IsCover nd.perm a b then some ⟨rswap nd.perm a b, (nd.perm b : ℕ)⟩ else none
    else []

/-- The operations performed by `stepNode`, one per element visited by each of its three
nested traversals: the scan building `A`, the scan over `b`, and, for each surviving pair,
the value comparison plus the scan between `a` and `b` made by the cover test. -/
def stepCost (k : ℕ) (nd : Node n) : ℕ :=
  (scanA n k).length +
    ((scanB n k).map fun b =>
      1 + (if (nd.perm b : ℕ) < nd.bound then
             ((eligible k nd).map fun a => 1 + (between a b).length).sum
           else 0)).sum

/-- The frontier `F` after `t` layers, starting from `(u, ∞)`. -/
def frontier (k : ℕ) (u : Perm n) : ℕ → List (Node n)
  | 0 => [⟨u, n + 1⟩]
  | t + 1 => (frontier k u t).flatMap (stepNode k)

/-- The list `T` returned by `pieri u p k`: one pair `(w, t)` per admissible `k`-chain from
`u` of length `t ≤ p`, recording its endpoint. -/
def table (u : Perm n) (p k : ℕ) : List (Perm n × ℕ) :=
  (List.range (p + 1)).flatMap fun t => (frontier k u t).map fun nd => (nd.perm, t)

/-- `M := |T|`. -/
def M (u : Perm n) (p k : ℕ) : ℕ := (table u p k).length

/-- The total cost: the algorithm expands the frontiers of layers `0, …, p-1`. -/
def cost (u : Perm n) (p k : ℕ) : ℕ :=
  ((List.range p).map fun t => ((frontier k u t).map (stepCost k)).sum).sum

/-! ### What a step produces -/

theorem mem_stepNode {k : ℕ} {nd nd' : Node n} (h : nd' ∈ stepNode k nd) :
    ∃ a b : Fin n, (a : ℕ) < k ∧ k ≤ (b : ℕ) ∧ IsCover nd.perm a b ∧
      (nd.perm b : ℕ) < nd.bound ∧
      nd' = ⟨rswap nd.perm a b, (nd.perm b : ℕ)⟩ := by
  simp only [stepNode, List.mem_flatMap] at h
  obtain ⟨b, hb, hmem⟩ := h
  by_cases hbb : (nd.perm b : ℕ) < nd.bound
  · rw [if_pos hbb] at hmem
    simp only [List.mem_filterMap] at hmem
    obtain ⟨a, ha, hres⟩ := hmem
    by_cases hcov : IsCover nd.perm a b
    · rw [if_pos hcov] at hres
      exact ⟨a, b, mem_scanA (mem_eligible ha), mem_scanB hb, hcov, hbb,
        (Option.some.inj hres).symm⟩
    · rw [if_neg hcov] at hres
      exact absurd hres (by simp)
  · rw [if_neg hbb] at hmem
    exact absurd hmem (by simp)

/-- Each step strictly lowers the bound: `β₀ > β₁ > ⋯`.  This is the mechanism behind
Lemma `lemma:pierimultfree`. -/
theorem bound_lt_of_mem_stepNode {k : ℕ} {nd nd' : Node n} (h : nd' ∈ stepNode k nd) :
    nd'.bound < nd.bound := by
  obtain ⟨a, b, -, -, -, hbb, rfl⟩ := mem_stepNode h
  exact hbb

/-- Each step raises the length by exactly one and moves up in Bruhat order. -/
theorem step_len {k : ℕ} {nd nd' : Node n} (h : nd' ∈ stepNode k nd) :
    len nd'.perm = len nd.perm + 1 ∧ BruhatLE nd.perm nd'.perm := by
  obtain ⟨a, b, hak, hkb, hcov, -, rfl⟩ := mem_stepNode h
  have hab : a < b := Fin.lt_def.mpr (lt_of_lt_of_le hak hkb)
  have hlen := len_rswap nd.perm hab hcov
  exact ⟨hlen, BruhatLE.step (BruhatLE.refl _) hlen⟩

/-- Every node of the layer-`t` frontier has length `ℓ(u) + t` and lies above `u`. -/
theorem frontier_spec {k : ℕ} {u : Perm n} :
    ∀ (t : ℕ) {nd : Node n}, nd ∈ frontier k u t →
      len nd.perm = len u + t ∧ BruhatLE u nd.perm := by
  intro t
  induction t with
  | zero =>
    intro nd h
    simp only [frontier, List.mem_singleton] at h
    subst h
    exact ⟨rfl, BruhatLE.refl _⟩
  | succ t ih =>
    intro nd h
    simp only [frontier, List.mem_flatMap] at h
    obtain ⟨nd₀, hnd₀, hstep⟩ := h
    obtain ⟨hlen₀, hbr₀⟩ := ih hnd₀
    obtain ⟨hlen, hbr⟩ := step_len hstep
    exact ⟨by rw [hlen, hlen₀]; ring, bruhat_trans hbr₀ hbr⟩

/-! ### The support bound -/

open scoped Classical in
/-- The set `{w ∈ Sₙ | w ≥ u, ℓ(w) ≤ ℓ(u) + p}` named by the proposition. -/
noncomputable def supportSet (u : Perm n) (p : ℕ) : Finset (Perm n) :=
  Finset.univ.filter fun w => BruhatLE u w ∧ len w ≤ len u + p

theorem mem_supportSet {u w : Perm n} {p : ℕ} :
    w ∈ supportSet u p ↔ BruhatLE u w ∧ len w ≤ len u + p := by
  classical
  simp [supportSet]

/-- Every endpoint returned by `pieri u p k` lies above `u` in Bruhat order and has length
at most `ℓ(u) + p`. -/
theorem endpoints_mem_supportSet (u : Perm n) (p k : ℕ) :
    ∀ w ∈ (table u p k).map Prod.fst, w ∈ supportSet u p := by
  intro w hw
  simp only [table, List.mem_map, List.mem_flatMap, List.mem_range] at hw
  obtain ⟨q, ⟨t, ht, nd, hnd, rfl⟩, rfl⟩ := hw
  obtain ⟨hlen, hbr⟩ := frontier_spec t hnd
  exact mem_supportSet.mpr ⟨hbr, by rw [hlen]; omega⟩

/-- **Multiplicity-freeness** (Lemma `lemma:pierimultfree` of `writing/schubmult.tex`).
Distinct admissible `k`-chains from `u` have distinct endpoints, so the endpoint list
returned by `pieri` is duplicate-free.

For chains of a common length this is the multiplicity-freeness of the `e_t`-Pieri expansion,
[Bergeron–Sottile, *Skew Schubert functions and the Pieri formula for flag
manifolds*][bergeron2002skew]; endpoints of different lengths are distinct because a
`k`-strip cover raises Bruhat length by one, which is `step_len` above.

Assumed, not proved: it is the one piece of established combinatorics this file takes from
the literature, and it is what `#print axioms` reports.  The `#guard`s in the `Guards`
section check it on 99 instances covering 562 chains. -/
axiom table_nodup (u : Perm n) (p k : ℕ) : ((table u p k).map Prod.fst).Nodup

/-- A duplicate-free list drawn from a finite set is no longer than that set. -/
theorem length_le_card_of_nodup {α : Type*} {l : List α} {S : Finset α}
    (hnd : l.Nodup) (hsub : ∀ x ∈ l, x ∈ S) : l.length ≤ S.card := by
  classical
  rw [← List.toFinset_card_of_nodup hnd]
  exact Finset.card_le_card fun x hx => hsub x (List.mem_toFinset.mp hx)

/-- **Proposition `prop:piericomplexity`, support half.**
`M ≤ #{w ∈ Sₙ | w ≥ u, ℓ(w) ≤ ℓ(u) + p}`. -/
theorem M_le_card (u : Perm n) (p k : ℕ) : M u p k ≤ (supportSet u p).card := by
  have hlen : (table u p k).length = ((table u p k).map Prod.fst).length := by simp
  rw [M, hlen]
  exact length_le_card_of_nodup (table_nodup u p k) (endpoints_mem_supportSet u p k)

/-! ### The cost bound -/

private theorem list_sum_le {l : List ℕ} {c : ℕ} (h : ∀ x ∈ l, x ≤ c) :
    l.sum ≤ l.length * c := by
  induction l with
  | nil => simp
  | cons x xs ih =>
    have h1 : x ≤ c := h x (by simp)
    have h2 : xs.sum ≤ xs.length * c := ih fun y hy => h y (by simp [hy])
    simp only [List.sum_cons, List.length_cons]
    calc x + xs.sum ≤ c + xs.length * c := Nat.add_le_add h1 h2
      _ = (xs.length + 1) * c := by ring

private theorem sum_map_mul (l : List ℕ) (f : ℕ → ℕ) (c : ℕ) :
    (l.map fun i => f i * c).sum = (l.map f).sum * c := by
  induction l with
  | nil => simp
  | cons x xs ih => simp [ih, Nat.add_mul]

/-- Cost of processing one frontier node, in the paper's shape: `O(k)` to build `A`, then
`n - k` choices of `b` against at most `k` elements of `A`, each pair costing `O(n)`. -/
def nodeCost (n k : ℕ) : ℕ := k + (n - k) * (1 + k * (n + 1))

/-- **The per-node cost bound.**  Each summand of `stepCost` is bounded by the length of the
list `stepNode` traverses to produce it. -/
theorem stepCost_le (k : ℕ) (nd : Node n) : stepCost k nd ≤ nodeCost n k := by
  have hinner : ∀ b : Fin n,
      (1 + (if (nd.perm b : ℕ) < nd.bound then
              ((eligible k nd).map fun a => 1 + (between a b).length).sum
            else 0)) ≤ 1 + k * (n + 1) := by
    intro b
    refine Nat.add_le_add_left ?_ 1
    split
    · refine le_trans (list_sum_le (c := n + 1) ?_) ?_
      · intro x hx
        simp only [List.mem_map] at hx
        obtain ⟨a, -, rfl⟩ := hx
        have := length_between_le a b
        omega
      · rw [List.length_map]
        exact Nat.mul_le_mul_right _
          ((length_eligible_le k nd).trans (length_scanA_le k))
    · exact Nat.zero_le _
  refine Nat.add_le_add (length_scanA_le k) ?_
  refine le_trans (list_sum_le (c := 1 + k * (n + 1)) ?_) ?_
  · intro x hx
    simp only [List.mem_map] at hx
    obtain ⟨b, -, rfl⟩ := hx
    exact hinner b
  · rw [List.length_map]
    exact Nat.mul_le_mul_right _ (length_scanB_le k)

/-- The per-node cost is `O(n³)`. -/
theorem nodeCost_le_cube (k : ℕ) (hk : k ≤ n) : nodeCost n k ≤ n ^ 3 + n ^ 2 + 2 * n := by
  have h1 : n - k ≤ n := Nat.sub_le _ _
  have h2 : 1 + k * (n + 1) ≤ 1 + n * (n + 1) :=
    Nat.add_le_add_left (Nat.mul_le_mul_right _ hk) 1
  calc nodeCost n k = k + (n - k) * (1 + k * (n + 1)) := rfl
    _ ≤ n + n * (1 + n * (n + 1)) := Nat.add_le_add hk (Nat.mul_le_mul h1 h2)
    _ = n ^ 3 + n ^ 2 + 2 * n := by ring

/-- The frontiers expanded by the algorithm, summed over layers `0, …, p-1`. -/
theorem sum_frontier_le_M (u : Perm n) (p k : ℕ) :
    ((List.range p).map fun t => (frontier k u t).length).sum ≤ M u p k := by
  have hM : M u p k = ((List.range (p + 1)).map fun t => (frontier k u t).length).sum := by
    simp [M, table, List.length_flatMap]
  rw [hM, List.range_succ, List.map_append, List.sum_append]
  exact Nat.le_add_right _ _

/-- **Proposition `prop:piericomplexity`, cost half.**  The algorithm runs in
`O(k (n-k) n · M)` elementary operations. -/
theorem cost_le_nodeCost_mul_M (u : Perm n) (p k : ℕ) :
    cost u p k ≤ nodeCost n k * M u p k := by
  have hlayer : ∀ t : ℕ, ((frontier k u t).map (stepCost k)).sum
      ≤ (frontier k u t).length * nodeCost n k := by
    intro t
    refine le_trans (list_sum_le (c := nodeCost n k) ?_) ?_
    · intro x hx
      simp only [List.mem_map] at hx
      obtain ⟨nd, -, rfl⟩ := hx
      exact stepCost_le k nd
    · rw [List.length_map]
  calc cost u p k
      = ((List.range p).map fun t => ((frontier k u t).map (stepCost k)).sum).sum := rfl
    _ ≤ ((List.range p).map fun t => (frontier k u t).length * nodeCost n k).sum :=
        List.sum_le_sum fun t _ => hlayer t
    _ = ((List.range p).map fun t => (frontier k u t).length).sum * nodeCost n k :=
        sum_map_mul _ _ _
    _ ≤ M u p k * nodeCost n k := Nat.mul_le_mul_right _ (sum_frontier_le_M u p k)
    _ = nodeCost n k * M u p k := Nat.mul_comm _ _

/-- The cost bound in the `O(n³)`-per-enumerated-permutation form. -/
theorem cost_le_cube (u : Perm n) (p k : ℕ) (hk : k ≤ n) :
    cost u p k ≤ (n ^ 3 + n ^ 2 + 2 * n) * M u p k :=
  le_trans (cost_le_nodeCost_mul_M u p k)
    (Nat.mul_le_mul_right _ (nodeCost_le_cube k hk))

/-- **Proposition `prop:piericomplexity`.**  The cost of `pieri u p k` is `O(n³)` per
element of `{w ∈ Sₙ | w ≥ u, ℓ(w) ≤ ℓ(u) + p}`. -/
theorem pieri_complexity (u : Perm n) (p k : ℕ) (hk : k ≤ n) :
    cost u p k ≤ (n ^ 3 + n ^ 2 + 2 * n) * (supportSet u p).card :=
  le_trans (cost_le_cube u p k hk) (Nat.mul_le_mul_left _ (M_le_card u p k))

/-! ### Agreement with the reference implementation

The counts below were checked against `elem_sym_perms` in
`src/schubmult/utils/schub_lib.py`, which enumerates the same chains.  One-line notation is
`0`-indexed here and `1`-indexed there. -/

section Guards

set_option linter.hashCommand false

private def sw (a b : Fin 6) : Perm 6 := Equiv.swap a b

-- `u = 123456`
#guard M (1 : Perm 6) 2 2 = 3
#guard M (1 : Perm 6) 3 3 = 4
-- `u = 321456`
#guard M (sw 0 2) 2 3 = 7
-- `u = 143256`
#guard M (sw 1 3) 1 1 = 4
#guard M (sw 1 3) 2 2 = 6
-- `u = 241356`
#guard M (sw 0 1 * sw 2 3 * sw 1 2) 3 3 = 8

private def guardSamples : List (Perm 6) :=
  [ 1, sw 0 1, sw 0 2, sw 1 3, sw 2 5,
    sw 0 1 * sw 1 2, sw 0 3 * sw 1 2, sw 0 2 * sw 1 3,
    sw 0 1 * sw 2 3 * sw 1 2, sw 0 5 * sw 1 4 * sw 2 3,
    sw 0 1 * sw 1 2 * sw 2 3 * sw 3 4 ]

private def guardParams : List (ℕ × ℕ) :=
  [(1,1),(1,2),(2,2),(2,3),(3,3),(3,4),(4,4),(4,5),(5,5)]

-- `table_nodup` on a spread of instances: 99 cases, 562 chains.
#guard guardSamples.all fun u =>
  guardParams.all fun pk => decide (((table u pk.1 pk.2).map Prod.fst).Nodup)

#guard (guardSamples.map fun u => (guardParams.map fun pk => M u pk.1 pk.2).sum).sum = 562

-- `len` is the inversion count.
#guard len (1 : Perm 6) = 0
#guard len (sw 0 1) = 1
#guard len (sw 1 3) = 3
#guard len (sw 0 5 * sw 1 4 * sw 2 3) = 15

end Guards

end Schubmult.Pieri
