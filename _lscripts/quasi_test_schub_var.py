from schubmult import *

def varpull_schub(perm, varnum, tvar):
    ring = DSx([], "t").ring
    result = ring.zero
    for new_perm, L in pull_out_var(varnum, perm):
        coeff = prod([tvar[varnum] - tvar[i] for i in L])
        if coeff != 0:
            result += coeff * ring(new_perm)
    return result

#really a divdiff
def check_trim_gpos():
    pass