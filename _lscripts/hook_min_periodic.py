"""Check that log F(n) - n log gamma + log(n+1) is asymptotically periodic in log2 n."""

import math

import numpy as np

N = 30000
lF = np.zeros(N + 1)
for n in range(1, N + 1):
    m = (n + 1) // 2
    lF[n] = math.log(n) + np.min(lF[0:m] + lF[n - m:n][::-1])

lgamma = sum(math.log(2 ** m - 1) / 2 ** m for m in range(1, 200))
res = lF - np.arange(N + 1) * lgamma + np.log(np.arange(N + 1) + 1.0)

print("dyadic window    min residual    max residual    amplitude")
k = 5
while 2 ** (k + 1) <= N:
    lo, hi = 2 ** k, min(2 ** (k + 1), N)
    r = res[lo:hi + 1]
    print(f"[2^{k:2d}, 2^{k + 1:2d}]    {r.min():11.7f}   {r.max():11.7f}   {r.max() - r.min():9.7f}")
    k += 1

print(f"\nlog(gamma/4) = {lgamma - math.log(4):.7f}   (predicted residual at n = 2^h - 1)")
print("  n=2^h-1: ", [f"{res[2 ** h - 1]:.7f}" for h in range(10, 15)])
print("  n=2^h:   ", [f"{res[2 ** h]:.7f}" for h in range(10, 15)])

for k in [11, 12, 13, 14]:
    lo, hi = 2 ** k, 2 ** (k + 1)
    b = lo + int(np.argmax(res[lo:hi + 1]))
    print(f"argmax in [2^{k},2^{k+1}] at n={b}  n/2^k={b / 2 ** k:.5f}  res={res[b]:.7f}")

print(f"\nC(n)=exp(res) in [{math.exp(res[1024:].min()):.5f}, {math.exp(res[1024:].max()):.5f}] for n >= 1024")
print(f"gamma   = {math.exp(lgamma):.10f}")
print(f"e*gamma = {math.e * math.exp(lgamma):.10f}")
