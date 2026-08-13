import math
from code_product_gap import Lam_from_code

print("rectangular codes lam = (c^k) in S_{k+c}:  Lambda vs binom(k+c,c) vs box (c+1)^k")
for c in range(1, 5):
    row = []
    for k in range(1, 7):
        N = k + c
        L = Lam_from_code(tuple([c] * k), N)
        b = math.comb(k + c, c)
        box = (c + 1) ** k
        row.append(f"k={k}: L={L} C={b} {'OK' if L == b else 'NO'} box={box} r={box / L:.1f}")
    print(f" c={c}: " + " | ".join(row))

print()
print("all-ones code (1^k) in S_{k+1}: Lambda=k+1, box=2^k")
for k in range(1, 12):
    L = Lam_from_code(tuple([1] * k), k + 1)
    print(f"  k={k:2d}  Lambda={L:3d}  box={2 ** k:5d}  ratio={2 ** k / L:8.2f}")
