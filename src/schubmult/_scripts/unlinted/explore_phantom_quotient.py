"""Explore the phantom quotient algebra.

Usage:
    conda activate schubmult_312
    python src/schubmult/_scripts/unlinted/explore_phantom_quotient.py [size]

Shows:
  1. Which permutations have phantom CEM terms at given size
  2. How phantoms identify RC graphs from different permutations
  3. Whether phantom contributions cancel in products
"""
import sys
import itertools

from schubmult.combinatorics.permutation import Permutation
from schubmult.combinatorics.rc_graph import RCGraph
from schubmult.rings.combinatorial.phantom_quotient_algebra import PhantomQuotientAlgebra

size = int(sys.argv[1]) if len(sys.argv) > 1 else 3

q = PhantomQuotientAlgebra(size)
perms = [p for p in Permutation.all_permutations(size + 1) if p.inv > 0]

print(f"=== Phantom Quotient Algebra at size {size} ===\n")

# 1. Which perms have phantoms?
print("--- Permutations with phantom CEM terms ---")
phantom_perms = []
for perm in perms:
    gen = q.phantom_generators(perm)
    if len(gen) > 0:
        phantom_perms.append(perm)
        print(f"  {perm} (code={perm.trimcode}): {len(gen)} phantom tensor terms")

if not phantom_perms:
    print("  None found at this size.\n")
    # Try larger size
    print(f"  Trying size {size + 1}...")
    q2 = PhantomQuotientAlgebra(size + 1)
    perms2 = [p for p in Permutation.all_permutations(size + 2) if p.inv > 0]
    for perm in perms2:
        gen = q2.phantom_generators(perm)
        if len(gen) > 0:
            phantom_perms.append(perm)
            print(f"  {perm} (code={perm.trimcode}): {len(gen)} phantom tensor terms")
    if not phantom_perms:
        print("  Still none. Phantoms appear at larger sizes.")
        sys.exit(0)
    q = q2
    size = size + 1
    perms = perms2

print()

# 2. Identification report for first few phantom perms
print("--- Phantom identification details ---")
for perm in phantom_perms[:3]:
    q.phantom_report(perm)
    print()

# 3. Schub elements: full vs non-phantom
print("--- Schubert elements: full vs non-phantom ---")
for perm in phantom_perms[:3]:
    s_full = q.schub_elem(perm)
    s_clean = q.schub_elem_nonphantom(perm)
    diff = q.from_dict({k: s_full.get(k, 0) - s_clean.get(k, 0) for k in set(s_full) | set(s_clean)})
    print(f"  S_{perm}:")
    print(f"    Full:  {dict(s_full)}")
    print(f"    Clean: {dict(s_clean)}")
    if any(v != 0 for v in diff.values()):
        print(f"    Diff:  {dict(diff)}  <-- phantoms shifted RC graph coefficients!")
    else:
        print(f"    Diff:  zero (phantoms cancel in Schubert element itself)")
    print()

# 4. Products: do phantoms cancel?
print("--- Product analysis ---")
tested = 0
for perm1, perm2 in itertools.product(phantom_perms[:3], repeat=2):
    full, clean, diff = q.ideal_intersection_report(perm1, perm2)
    tested += 1
    print()

if tested == 0:
    print("  No phantom perms to test products with.")

# 5. Ring structure: what elements generate the phantom ideal?
print("--- Phantom ideal generators ---")
for perm in phantom_perms[:3]:
    gen = q.phantom_generators(perm)
    print(f"  Phantom generator for {perm}:")
    for key, coeff in gen.items():
        rc = q._brc.key_to_rc_graph(key)
        print(f"    coeff={coeff}, squashes_to perm={rc.perm}, key_factors={len(key)}")
    print()

print("Done.")
