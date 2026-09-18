from schubmult import *
from schubmult.rings.polynomial_algebra import *
from schubmult.symbolic import sympify, S
from schubmult.symbolic.poly.variables import GeneratingSet, CustomGeneratingSet
from schubmult.rings.schubert.double_schubert_ring import DoubleSchubertRing
from schubmult.combinatorics.double_forest import *

t = GeneratingSet("t")


def varshiftsub(pol, arr, index, subvar):
    shifted = pol.subs(arr[index], subvar).subs({arr[i]: arr[i - 1] for i in range(index + 1, len(arr) - 1)})
    return shifted

def actual_quasi_dd(pol, arr, i, tvar):
    from sympy import quo
    workpol = quo((varshiftsub(pol, arr, i + 1, tvar) - varshiftsub(pol, arr, i, tvar)),arr[i] - tvar)
    return workpol

def schub_quasi_dd(perm, arr, i, tvar):
    cgs = CustomGeneratingSet([*arr[:i+1], tvar, *arr[i+1:]])
    ring = DoubleSchubertRing(cgs, t)
    if perm[i - 1] < perm[i]:
        return 0
    the_dict = ring(perm.swap(i - 1, i)).pull_out_gen(tvar)
    return the_dict


def schub_quasi_dd_varset(perm, arr, i, tvar):
    cgs = CustomGeneratingSet([*arr[:i+1], tvar, *arr[i+1:]])
    ring = DoubleSchubertRing(cgs, t)
    if perm[i - 1] < perm[i]:
        return 0
    the_dict = ring(perm.swap(i - 1, i))
    return the_dict, ring.genset


def full_schub_quasi_dd(elem, arr, i, tvar):
    result = S.Zero
    if elem == 0:
        return S.Zero
    for perm, coeff in elem.items():
        result += coeff * schub_quasi_dd(perm, arr, i, tvar)
    return result



def subset_from_seq(seq):
    if len(seq) <= 1:
        return tuple(seq)
    for i in range(len(seq) - 1):
        if seq[i] >= seq[i + 1]:
            seq2 = [*seq[:i], seq[i+1], seq[i] + 1, *seq[i + 2:]]
            return subset_from_seq(seq2)
    return tuple(seq)

def vartrim(pol, arr, index):
    return pol.subs({arr[i]: arr[i + 1] for i in range(index, len(arr) - 1)})


def t_index_iA(i, A):
    """Return j such that t_{i, A} = t_j, i.e. j = (\\bar A)_i (1-indexed i-th elt of N \\ A)."""
    A_set = set(A)
    count = 0
    j = 0
    while True:
        j += 1
        if j not in A_set:
            count += 1
            if count == i:
                return j


def ev_A(pol, x_arr, A, t_gen, n_xvars):
    """ev_A f = f(t_{(bar A)_1}, ..., t_{(bar A)_{n_xvars}}; t) — substitute x_j -> t_{(bar A)_j}."""
    subs_dict = {x_arr[j]: t_gen[t_index_iA(j, A)] for j in range(1, n_xvars + 1)}
    return sympify(pol).subs(subs_dict)


def factorization_from_code(code):
    """Columns factorization of forest F=1^{c_1}·2^{c_2}·... in Thompson monoid."""
    seq = []
    for i, c in enumerate(code, start=1):
        seq.extend([i] * c)
    return seq


def a_F_polynomial(f_poly, x_arr, t_gen, code, n_xvars):
    """Compute a_F(t) = [ev ⋆ E_F] f_poly  for forest with given code via direct
    polynomial divided differences.

    For factorization (i_1, ..., i_k) of F:
        a_F = ev_{A_0} ∘ E_{i_1, A_1} ∘ ... ∘ E_{i_k, ∅}  f
    where A_j = i_{j+1} ⋆ ... ⋆ i_k.
    """
    from schubmult.mult.double import schubmult_double
    seq = factorization_from_code(code)
    cur = sympify(f_poly)
    tvars = []
    for j in range(len(seq) - 1, -1, -1):
        i_j = seq[j]
        A_j = subset_from_seq(seq[j + 1:])
        tvar = t_gen[t_index_iA(i_j, A_j)]
        #cur = actual_quasi_dd(cur, x_arr, i_j, tvar)
        tvars = [tvar] + tvars
    A_0 = subset_from_seq(seq)
    return ev_A(cur, x_arr, A_0, t_gen, n_xvars)


def schub_a_F_polynomial(schub, x_arr, t_gen, code, n_xvars):
    """Compute a_F(t) = [ev ⋆ E_F] f_poly  for forest with given code via direct
    polynomial divided differences.

    For factorization (i_1, ..., i_k) of F:
        a_F = ev_{A_0} ∘ E_{i_1, A_1} ∘ ... ∘ E_{i_k, ∅}  f
    where A_j = i_{j+1} ⋆ ... ⋆ i_k.
    """
    from schubmult.mult.double import schubmult_double
    seq = factorization_from_code(code)
    #cur = sympify(f_poly)
    cur  = schub
    #tvars = []
    for j in range(len(seq) - 1, -1, -1):
        i_j = seq[j]
        A_j = subset_from_seq(seq[j + 1:])
        tvar = t_gen[t_index_iA(i_j, A_j)]
        cur = full_schub_quasi_dd(cur, x_arr, i_j, tvar)
        # tvars = [*tvars, tvar]
    A_0 = subset_from_seq(seq)
    # the_perm = ~(uncode([a - 2 if a >= 2 else 0 for a in reversed(A_0)]))
    # the_get_perm = ~(uncode(tuple((a - 1 if a >= 1 else 0 for a in reversed(A_0)))))
    # result = 0
    # for u, coeff in schub.items():    
    #     coeff2 = schubmult_double({the_perm: 1}, u, [*tvars, *x_arr[1:]], t_gen).get(the_get_perm, 0)
        # result += coeff * ev_A(coeff2, x_arr, A_0, t_gen, n_xvars)
    #print(coeff)
    # return result
    return ev_A(cur.as_polynomial() if isinstance(cur, DoubleSchubertElement) else sympify(cur), x_arr, A_0, t_gen, n_xvars)

def lrcoeff_a_F_polynomial(schub, x_arr, t_gen, code, n_xvars):
    """Compute a_F(t) = [ev ⋆ E_F] f_poly  for forest with given code via direct
    polynomial divided differences.

    For factorization (i_1, ..., i_k) of F:
        a_F = ev_{A_0} ∘ E_{i_1, A_1} ∘ ... ∘ E_{i_k, ∅}  f
    where A_j = i_{j+1} ⋆ ... ⋆ i_k.
    """
    from schubmult.mult.double import schubmult_double
    seq = subset_from_seq(factorization_from_code(code))
    #cur = sympify(f_poly)
    cur  = schub
    tvars = []
    for j in range(len(seq) - 1, -1, -1):
        i_j = seq[j]
        A_j = subset_from_seq(seq[j + 1:])
        tvar = t_gen[t_index_iA(i_j, A_j)]
        #cur = full_schub_quasi_dd(cur, x_arr, i_j, tvar)
        tvars = [tvar, *tvars]
    A_0 = subset_from_seq(seq)
    the_perm = ~(uncode([a - 1 if a >= 1 else 0 for a in reversed(A_0)]))
    the_get_perm = ~(uncode(tuple(reversed(A_0))))
    result = 0
    subs_dict = {x_arr[i + 1]: tvars[i] for i in range(len(tvars))}
    subs_dict2 = {x_arr[i]: 0 for i in range(len(tvars) + 1)} | {x_arr[len(tvars) + i + 1]: x_arr[i + 1] for i in range(len(x_arr) - len(tvars) - 1)}
    for u, coeff in schub.items():    
        coeff2 = schubmult_double({the_perm: 1}, u, x_arr, t_gen).get(the_get_perm, 0)
        result += coeff * ev_A(sympify(coeff2).subs(subs_dict).subs(subs_dict2), x_arr, A_0, t_gen, n_xvars)
    #print(coeff)
    return result
    #return ev_A(cur.as_polynomial() if isinstance(cur, DoubleSchubertElement) else sympify(cur), x_arr, A_0, t_gen, n_xvars)

def enum_forest_codes(length, max_sum):
    """All forest codes (c_1, ..., c_length) with sum <= max_sum, c_i >= 0."""
    if length == 0:
        yield ()
        return
    for c1 in range(max_sum + 1):
        for tail in enum_forest_codes(length - 1, max_sum - c1):
            yield (c1, *tail)



if __name__ == "__main__":
    from schubmult.abc import x
    import sys
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    DSt = DoubleSchubertRing(x, t)
    DForest = PolynomialAlgebra(DoubleForestPolyBasis(x, t))

    for perm in perms:
        if perm.inv == 0:
            continue
        schub = DSt(perm)
        # for i in perm.descents(zero_indexed=False):
        #     the_quasi_dd = schub_quasi_dd(perm, x, 1, t[1])
        #     real_quasi_dd = actual_quasi_dd(schub.expand(), x, 1, t[1])

        #     if the_quasi_dd == real_quasi_dd or the_quasi_dd.almosteq(real_quasi_dd):
        #         print(f"Success for {perm} at descent {i}")
        #     else:
        #         print(f"Failure for {perm} at descent {i}")
        #         print(f"Computed: {the_quasi_dd}")
        #         print(f"Expected: {real_quasi_dd}")
        # ⋆-monoid forest-coefficient extraction (Section 10, arXiv:2504.15234):
        #
        #   S_w = sum_F a_F(t) P_F   with   a_F(t) = [ev ⋆ E_{i_1} ⋆ ... ⋆ E_{i_k}] S_w
        #
        # for ANY Thompson factorization (i_1, ..., i_k) of F. The right-hand side equals
        #
        #   ev_{A_0} ∘ E_{i_1, A_1} ∘ ... ∘ E_{i_{k-1}, A_{k-1}} ∘ E_{i_k, ∅}  S_w
        #
        # with A_j = i_{j+1} ⋆ ... ⋆ i_k, where the ⋆-monoid relation i ⋆ j = j ⋆ (i+1)
        # for i ≥ j identifies a sequence with its canonical strictly-increasing form
        # (interpreted as a finite subset A ⊂ N), and t_{i, A} := t_{(\bar A)_i} where
        # \bar A = N \ A. The bug in the previous attempt was indexing t into A directly
        # rather than into its complement \bar A.
        #
        # We enumerate forests via codes of length n-1 with sum ≤ inv(w) (since a_F has
        # t-degree inv(w) - |F| ≥ 0), use the columns factorization 1^{c_1}·2^{c_2}·...
        # for each, and compute a_F by direct polynomial divided differences. This avoids
        # the over-counting of the dict/BFS approach (different factorizations of F all
        # give the same a_F so summing over factorizations multiplies by #Trim(F)).
        
        forest_dict = sympify(0)
        for cd in enum_forest_codes(n - 1, perm.inv):
            a = schub_a_F_polynomial(schub, x, t, cd, n)
            #a = lrcoeff_a_F_polynomial(schub, x, t, cd, n)
            a = sympify(a)
            if a != 0:
                forest_dict += a * DForest(*cd)

        if forest_dict == 0:
            result = DForest.zero
        else:
            result = forest_dict
            # forest_dict is already a sum of a_F·DForest(cd) terms; coerce to DForest element
            result = forest_dict
        result2 = DForest.from_expr(schub.expand(), length=n - 1)
        # diff_expr = (sympify(result.expand() if hasattr(result, "expand") else result)
        #              - sympify(result2.expand())).expand()
        if result != result2 and not result.almosteq(result2):
            print(f"Forest FAIL for {perm}")
            print(f"  got:      {result}")
            print(f"  expected: {result2}")
            #print(f"  diff:     {diff_expr}")
        else:
            print(f"Forest OK   for {perm}")
            print(f"  got:      {result}")
        