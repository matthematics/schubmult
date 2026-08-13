import math
import numpy as np

N = 30000
lF = np.zeros(N + 1)
for n in range(1, N + 1):
    m = (n + 1) // 2
    lF[n] = math.log(n) + np.min(lF[0:m] + lF[n - m:n][::-1])

lgamma = sum(math.log(2 ** m - 1) / 2 ** m for m in range(1, 200))
idx = np.arange(N + 1)
res = lF - idx * lgamma + np.log(idx + 1.0)
c0 = lgamma - math.log(4)

print(f"log(gamma/4) = {c0:.10f}")
print(f"min residual over 1..{N}: {res[1:].min():.10f}  at n={int(np.argmin(res[1:])) + 1}")
print(f"all residuals >= log(gamma/4)? {bool((res[1:] >= c0 - 1e-12).all())}")
print(f"max residual over 1..{N}: {res[1:].max():.10f}  at n={int(np.argmax(res[1:])) + 1}")

# growth rate of the max residual
xs, ys = [], []
k = 6
while 2 ** (k + 1) <= N:
    lo, hi = 2 ** k, min(2 ** (k + 1), N)
    xs.append(math.log(2 ** k))
    ys.append(res[lo:hi + 1].max())
    k += 1
A = np.polyfit(xs, ys, 1)
print(f"\nmax residual per window ~ {A[0]:.5f} * ln n + {A[1]:.5f}")
print(f"   => F(n) = gamma^(n+1)/(4(n+1)) * n^delta,  delta in [0, {A[0]:.4f}]")

print(f"\ngamma      = {math.exp(lgamma):.12f}")
print(f"e*gamma    = {math.e * math.exp(lgamma):.12f}")
print(f"log2 gamma = {lgamma / math.log(2):.12f}")

# exact max linear extension counts, small n
Fex = [1] * 21
for n in range(1, 21):
    Fex[n] = n * min(Fex[i] * Fex[n - 1 - i] for i in range(n))
print("\n n   F(n)        max_T e(T) = n!/F(n)")
for n in range(1, 21):
    print(f"{n:2d}  {Fex[n]:10d}   {math.factorial(n) // Fex[n]}")
