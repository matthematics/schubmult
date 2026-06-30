# Transplant Map (Wise First Pass)

This map migrates content from writing/schubert_bialgebra.tex into writing/schubert_bialgebra_paperA_mockup.tex in a low-risk order.

Principle: move proof-critical material first, then only add upstream dependencies if compilation or logic requires them.

## Target Anchors (Mockup)

- Introduction: line 44
- Main theorem: line 53
- Minimal Algebraic Setup: line 77
- Bounded RC graphs + row-factorization theorem: lines 92-107
- Forest Quotient Layer: line 109
- Lift Product and Counting Mechanism: line 127
- Proof of Main Theorem: line 151
- Geometric Corollary: line 162

## Phase 0: Lock Scope (No Copy Yet)

- Keep this pass focused on forest LR theorem and immediate machinery.
- Do not import dual key/dual slide sections unless forced by a missing dependency.
- Do not import notation catalog or appendix proofs in this pass.

Stop check:
- Confirm intro claims match reduced scope before any detailed proof import.

## Phase 1: Import Main Statement + Endgame Core

Target:
- writing/schubert_bialgebra_paperA_mockup.tex (Introduction and Proof section)

Source blocks (primary):
- Main forest theorem statement: writing/schubert_bialgebra.tex:180
- Forest LR endgame subsection: writing/schubert_bialgebra.tex:3839
- Witness/technical chain in endgame:
  - proposition:liftshufflewords at writing/schubert_bialgebra.tex:3848
  - lemma:foreststable at writing/schubert_bialgebra.tex:3874
  - lemma:forestlead at writing/schubert_bialgebra.tex:3901
  - lemma:forestchoice at writing/schubert_bialgebra.tex:3908
  - dual forest coproduct corollary (optional in pass 1): writing/schubert_bialgebra.tex:3988

How to transplant:
- Replace placeholder theorem at mockup line 53 with near-final statement from source line 180.
- Replace placeholder proposition/proof at mockup lines 154-160 with real extraction argument from source lines 3839+.

Stop check:
- Proof of main theorem in mockup must be self-contained modulo explicit forward references only.

## Phase 2: Import Lift Product Spine

Target:
- Lift Product and Counting Mechanism (mockup line 127 onward)

Source blocks:
- Multiplactic setup: writing/schubert_bialgebra.tex:3405
- Critical statements:
  - proposition:boxproduct at writing/schubert_bialgebra.tex:3417
  - theorem:placticiso at writing/schubert_bialgebra.tex:3451
  - theorem:schubmoduleiso at writing/schubert_bialgebra.tex:3464
  - proposition:eproductworks at writing/schubert_bialgebra.tex:3495
  - theorem:schubproduct at writing/schubert_bialgebra.tex:3514
- Lift product section:
  - section start writing/schubert_bialgebra.tex:3532
  - proposition:elemfactor at writing/schubert_bialgebra.tex:3536
  - lemma:liftweight at writing/schubert_bialgebra.tex:3571
  - lemma:liftcommute at writing/schubert_bialgebra.tex:3610
  - theorem at writing/schubert_bialgebra.tex:3621
- Slide quotient subset only if needed for forest proof chain:
  - proposition:quasicrystalaxioms at writing/schubert_bialgebra.tex:3691
  - theorem:quasicharacter at writing/schubert_bialgebra.tex:3728
  - lemma:quasitensorswap at writing/schubert_bialgebra.tex:3739
  - lemma:quasiprod_plactic at writing/schubert_bialgebra.tex:3746
  - lemma:liftequivariant at writing/schubert_bialgebra.tex:3763
  - proposition:quasiprod at writing/schubert_bialgebra.tex:3774
  - theorem:slideproduct at writing/schubert_bialgebra.tex:3826

How to transplant:
- First import only the minimum statements and short proof skeletons referenced by Phase 1.
- Delay long quasi-crystal discussions until a missing link is identified.

Stop check:
- Every forward reference in forest endgame must now resolve to a local statement in Paper A.

## Phase 3: Import Forest Quotient Layer

Target:
- Forest Quotient Layer (mockup line 109 onward)

Source blocks:
- Forest classes section start: writing/schubert_bialgebra.tex:3020
- Required items:
  - proposition:forestdualcoefficients at writing/schubert_bialgebra.tex:3102
  - theorem:forestformula at writing/schubert_bialgebra.tex:3132
  - definition:gamma at writing/schubert_bialgebra.tex:3160
  - theorem:LRforest at writing/schubert_bialgebra.tex:3198

How to transplant:
- Keep only definitions/results directly used by Phase 1-2 arguments.
- Move examples and figures in this section to a deferred appendix note.

Stop check:
- The quotient compatibility claim in mockup line 147 is now fully justified.

## Phase 4: Import Minimal RC-Ring Infrastructure

Target:
- Minimal Algebraic Setup (mockup line 77 onward)

Source blocks:
- RC graph and bounded module:
  - section:brc at writing/schubert_bialgebra.tex:662
  - definition:rcgraph at writing/schubert_bialgebra.tex:666
- Zero/clip/trim infrastructure:
  - definition:zeromap at writing/schubert_bialgebra.tex:1065
  - theorem:zero_bijection at writing/schubert_bialgebra.tex:1323
  - section:ringproduct at writing/schubert_bialgebra.tex:1397
- Use proposition:pullindex only if explicitly required by imported proof text:
  - writing/schubert_bialgebra.tex:350 and appendix proof writing/schubert_bialgebra.tex:4184

How to transplant:
- Keep this section terse: definitions + one structural theorem statement + references.
- Avoid importing long algorithm prose unless proof steps explicitly depend on it.

Stop check:
- Row-factorization theorem at mockup line 101 has a complete proof path in Paper A.

## Phase 5: Geometric Corollary and Intro Polish

Target:
- Geometric Corollary section (mockup line 162) and introduction (line 44)

Source blocks:
- Geometric interpretation subsection: writing/schubert_bialgebra.tex:4001
- Cup product corollary: writing/schubert_bialgebra.tex:4011

How to transplant:
- Keep one corollary + short proof only.
- Move branching-comorphism discussion to a short remark or companion note.

Stop check:
- No intro claim references results omitted from Paper A.

## Defer to Companion Paper (First-Pass Move-Out)

- RC graph crystals dual-key package: writing/schubert_bialgebra.tex:2374-2558.
- Dual slide package and general LR template, unless needed by a retained proof:
  - writing/schubert_bialgebra.tex:3257-3377.
- Notation compendium and standalone appendix proofs not needed by retained chain:
  - writing/schubert_bialgebra.tex:4047 onward.

## Risk Notes (Why This Is Wise)

- The dependency list in writing/schubert_bialgebra_cut_plan.md is intentionally conservative and may over-include due to narrative references in the introduction.
- This phased order minimizes rework by importing theorem-endgame first, then only the prerequisites it truly needs.
- If a dependency fails, add one missing upstream lemma rather than bulk-importing an entire section.

## Mechanical Workflow Per Phase

1. Copy only labeled statements first.
2. Recompile.
3. Resolve undefined references.
4. Copy proof text only for unresolved statements.
5. Recompile again.

Repeat until phase stop check passes.
