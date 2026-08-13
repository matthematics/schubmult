"""Tight asymptotics for max # linear extensions of a binary search tree on n nodes.

e(T) = n! / prod_x h_T(x)   (Knuth's hook length formula for forests)
max_T e(T) = n! / F(n),  F(n) = min_T prod h = n * min_{i+j=n-1} F(i) F(j).
"""

import math

N = 4000

# exact small values
Fex = [1] * 25
for n in range(1, 25):
    Fex[n] = n * min(Fex[i] * Fex[n - 1 - i] for i in range(n))
print("exact F(n), n=0..24:")
print(Fex)

# log-domain DP
lF = [0.0] * (N + 1)
argmin = [0] * (N + 1)
for n in range(1, N + 1):
    best = None
    bi = 0
    for i in range((n + 1) // 2):  # symmetric
        v = lF[i] + lF[n - 1 - i]
        if best is None or v < best:
            best = v
            bi = i
    lF[n] = math.log(n) + best
    argmin[n] = bi

# gamma = prod_{m>=1} (2^m - 1)^(2^-m)
lgamma = sum(math.log(2 ** m - 1) / 2 ** m for m in range(1, 200))
gamma = math.exp(lgamma)
print(f"\ngamma = {gamma:.12f}   log2(gamma) = {lgamma / math.log(2):.12f}")

# perfect-tree closed form check
for h in range(1, 6):
    n = 2 ** h - 1
    p = 1
    for m in range(1, h + 1):
        p *= (2 ** m - 1) ** (2 ** (h - m))
    print(f"  h={h} n={n:5d}  perfect prod={p}  F(n)={Fex[n] if n < 25 else '-'}")

# is log F(n) - n log gamma + log n bounded?
print("\n   n   log2F/n    logF - n*lg + log(n+1)     balanced?")
for n in [7, 15, 31, 63, 100, 127, 200, 255, 500, 511, 1000, 1023, 2000, 2047, 4000]:
    d = lF[n] - n * lgamma + math.log(n + 1)
    print(f"{n:6d}  {lF[n] / n / math.log(2):8.5f}   {d:12.6f}")

# sweep the residual over a range to bound the oscillation
res = [lF[n] - n * lgamma + math.log(n + 1) for n in range(100, N + 1)]
print(f"\nresidual over n in [100,{N}]: min={min(res):.6f} max={max(res):.6f}")
print(f"  => F(n) = C(n) * gamma^n / (n+1), C(n) in [{math.exp(min(res)):.5f}, {math.exp(max(res)):.5f}]")

# where is the min attained: how balanced is the optimal split?
print("\nn, optimal left-subtree size i (n-1-i on right):")
for n in [10, 20, 50, 100, 200, 500, 1000, 2000, 4000]:
    print(f"  n={n:5d}  i={argmin[n]:5d}   (n-1)/2={(n - 1) / 2:.1f}")
