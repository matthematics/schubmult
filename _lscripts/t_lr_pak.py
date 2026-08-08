from schubmult import *
from schubmult.utils.tuple_utils import pad_tuple
from schubmult.mult.single import schubmult_py
from math import comb
import sympy

def t_value(t, lmd, mu, nu):
    if len(lmd) != len(nu) or len(lmd) != len(mu):
        raise ValueError("two partitions must have same length")
    if t == 0:
        return 1
    k = len(lmd)
    base_part_lm = [t * lmd[i] for i in range(k - 1, -1, -1)]
    base_part_mu = [t * mu[i] for i in range(k - 1, -1, -1)]
    base_part_nu = [t * nu[i] for i in range(k - 1, -1, -1)]
    base_part_b = [t * b for _ in range(k)]
    mu_perm = uncode(base_part_mu)
    right_mu_perm = mu_perm.minimal_dominant_above()
    nu_perm = uncode(base_part_nu)
    right_nu_perm = nu_perm * (~mu_perm) * right_mu_perm
    return schubmult_py({uncode(base_part_lm): 1}, right_mu_perm).get(right_nu_perm, 0)

#def interp_poly(t, values):


def diff_operator(n, fnc, val):
    if n == 0:
        return fnc(val)
    d1 = diff_operator(n - 1, fnc, val)
    d2 = diff_operator(n - 1, fnc, val + 1)
    return d2 - d1

def t_polynomial_coeff(lmd, nu, b, index):
    return diff_operator(index, lambda x: t_value(x, lmd, nu, b), 0)
    
    

def t_binomial_coeffs(lmd, mu, nu, zero_buffer=0, verify=True):
    """Coefficients c_k such that

        t_value(t, lmd, nu, b) = sum_k c_k * binomial(t, k),

    computed by finite differences (Newton forward formula):

        c_k = (Delta^k f)(0),   Delta f(x) = f(x + 1) - f(x).

    Since f is a polynomial of some degree d, Delta^k f == 0 for k > d.
    We keep computing until we see `zero_buffer` consecutive zeros, then
    (optionally) verify the reconstruction against fresh sample points.
    """
    f = lambda x: t_value(x, lmd, nu, b)
    coeffs = []
    k = 0
    zero_run = 0
    while zero_run <= zero_buffer:
        c = t_polynomial_coeff(lmd, mu, nu, k)
        coeffs.append(c)
        zero_run = zero_run + 1 if c == 0 else 0
        k += 1
    # drop trailing zeros -> coeffs has length deg(f) + 1
    while coeffs and coeffs[-1] == 0:
        coeffs.pop()

    if verify:
        # Newton forward reconstruction must reproduce f at fresh points.
        deg = len(coeffs) - 1
        for x in range(deg + 1, deg + 1 + zero_buffer + 2):
            recon = sum(coeffs[j] * comb(x, j) for j in range(len(coeffs)))
            if recon != f(x):
                raise ValueError(
                    f"binomial-basis reconstruction failed at t={x}: "
                    f"got {recon}, expected {f(x)} (try a larger zero_buffer)"
                )
    return coeffs


def t_binomial_poly(lmd, mu, nu, t=None, **kwargs):
    """Symbolic expansion sum_k c_k * binomial(t, k) as a sympy expression."""
    if t is None:
        t = sympy.Symbol("t")
    coeffs = t_binomial_coeffs(lmd, mu, nu, **kwargs)
    return sympy.Add(*[c * sympy.binomial(t, k) for k, c in enumerate(coeffs)])


def _partition(perm, k):
    return tuple(reversed(perm.code[:k]))


def main(n, k, b, t):
    import itertools
    comps = [tuple(reversed(perm.code[:k])) for perm in Permutation.all_permutations(n) if perm.descents() == {k - 1}]
    #bcomp = [b for _ in range(k)]
    for comp1, comp2 in itertools.product(comps, repeat=2):
        
        fishpot = {k for k, v in (Sx(uncode(comp1)) * Sx(uncode(comp2))).items() if v != 0}
        for pisser in fishpot:
            pispot = _partition(pisser, k)
            coeffs = t_binomial_coeffs(comp1, comp2, pispot)
            poly = t_binomial_poly(comp1, comp2, pispot, t)
            print(f"{comp1,comp2,pispot}: coeffs={coeffs}  ->  {poly}")
            #break 

if __name__ == "__main__":
    import sys
    t = sympy.Symbol("t")
    n = int(sys.argv[1])
    k = int(sys.argv[2])
    b = int(sys.argv[3])
    main(n, k, b, t)