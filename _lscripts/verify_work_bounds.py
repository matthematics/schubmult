"""Verify candidate factor bounds for the schubmult work W.

Measures, per (u,v):
  W        = exact inner-body executions (proven-faithful loop)
  L        = number of layers
  Phi_max  = max over layers of node count (up,vp pairs)
  B        = max over layers,up of PIERI branching |newperms|
  Delta    = max over layers,vp of kdown-triple count
  Vlab_max = max over layers of number of DISTINCT v-labels present
  distinct_up_max = max over layers of number of distinct up-perms

Tests the candidate upper bound  W <= L * Phi_max * B * Delta
and the node bound Phi_max <= distinct_up_max * Vlab_max,
and whether Vlab_max <= 2**l(v)  vs  2**l(mu).
"""

import itertools

from schubmult import Permutation, Sx
from schubmult.combinatorics.permutation import uncode
from schubmult.utils.schub_lib import compute_vpathdicts, elem_sym_perms


def measure(u, v):
    vn1 = ~v
    th = vn1.theta()
    if len(th) == 0 or th[0] == 0:
        return None
    mu = uncode(th)
    vmu = v * mu
    inv_vmu = vmu.inv
    inv_mu = mu.inv
    lmu = mu.inv
    th = list(th)
    while th and th[-1] == 0:
        th.pop()
    thL = len(th)
    vpd = compute_vpathdicts(th, vmu)
    mx_th = [0] * thL
    for index in range(thL):
        for vp in vpd[index]:
            for _v2, vdiff, _s in vpd[index][vp]:
                mx_th[index] = max(mx_th[index], th[index] - vdiff)

    inv_u = u.inv
    vpathsums = {u: {Permutation([1, 2]): 1}}
    W = 0
    Phi_max = 0
    B = 0
    Delta = 0
    Vlab_max = 0
    dup_max = 0
    for index in range(thL):
        newpathsums = {}
        N = sum(len(vpathsums[up]) for up in vpathsums)
        Phi_max = max(Phi_max, N)
        dup_max = max(dup_max, len(vpathsums))
        labels = set()
        for up in vpathsums:
            labels |= set(vpathsums[up].keys())
        Vlab_max = max(Vlab_max, len(labels))
        for up in vpathsums:
            inv_up = up.inv
            newperms = elem_sym_perms(up, min(mx_th[index], inv_mu - inv_vmu - (inv_up - inv_u)), th[index])
            B = max(B, len(newperms))
            for vp in vpathsums[up]:
                triples = vpd[index][vp]
                Delta = max(Delta, len(triples))
                for _v2, vdiff, _s in triples:
                    for _up2, udiff in newperms:
                        W += 1
                        if vdiff + udiff == th[index]:
                            newpathsums.setdefault(_up2, {})
                            newpathsums[_up2][_v2] = newpathsums[_up2].get(_v2, 0) + 1
        vpathsums = newpathsums

    return {"W": W, "L": thL, "Phi_max": Phi_max, "B": B, "Delta": Delta, "Vlab_max": Vlab_max, "dup_max": dup_max, "lv": v.inv, "lmu": lmu}


def main():
    import sys

    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    bound_ok = True
    node_ok = True
    vlab_lv_ok = True
    vlab_lmu_ok = True
    worst_slack = (1e9, None)
    checked = 0
    for n in range(2, nmax + 1):
        for u, v in itertools.product(Permutation.all_permutations(n), repeat=2):
            r = measure(u, v)
            if r is None:
                continue
            checked += 1
            bound = r["L"] * r["Phi_max"] * r["B"] * r["Delta"]
            if r["W"] > bound:
                bound_ok = False
                print("BOUND VIOLATED", list(u), list(v), r, "bound=", bound)
            else:
                if bound > 0 and r["W"] > 0:
                    slack = bound / r["W"]
                    if slack < worst_slack[0]:
                        worst_slack = (slack, (list(u), list(v), r, bound))
            if r["Phi_max"] > r["dup_max"] * r["Vlab_max"]:
                node_ok = False
                print("NODE BOUND VIOLATED", list(u), list(v), r)
            if r["Vlab_max"] > 2 ** r["lv"]:
                vlab_lv_ok = False
            if r["Vlab_max"] > 2 ** r["lmu"]:
                vlab_lmu_ok = False
    print(f"checked {checked} pairs up to n={nmax}")
    print(f"W <= L*Phi_max*B*Delta : {'HOLDS' if bound_ok else 'FAILS'}")
    print(f"Phi_max <= distinct_up * distinct_vlabels : {'HOLDS' if node_ok else 'FAILS'}")
    print(f"Vlab_max <= 2^l(v) : {'HOLDS' if vlab_lv_ok else 'FAILS'}")
    print(f"Vlab_max <= 2^l(mu) : {'HOLDS' if vlab_lmu_ok else 'FAILS'}")
    print(f"tightest bound/W slack: {worst_slack[0]:.2f} at {worst_slack[1][:2] if worst_slack[1] else None}")


if __name__ == "__main__":
    main()
