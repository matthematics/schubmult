if __name__ == "__main__":
    from schubmult import *
    import sys
    from sympy import pretty_print
    from schubmult.symbolic import S, expand
    from schubmult.abc import x
    from schubmult.utils.perm_utils import has_bruhat_descent

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    
    def dual_schub(perm, x):
        if perm.inv == 0:
            return S.One
        res = S.Zero
        for a in range(len(perm) - 1):
            for b in range(a + 1, len(perm)):
                if has_bruhat_descent(perm, a, b):
                    #res += ((x[a + 1]-x[a]) - (x[b + 1]-x[b])) * dual_schub(perm.swap(a, b), x)
                    res += sum([x[i] for i in range(a + 1, b + 1)]) * dual_schub(perm.swap(a, b), x)
        return res

    for perm in perms:
        ds = expand(dual_schub(perm, x))
        sch = ASx(perm, n).change_basis(WordBasis)
        print(perm.trimcode)
        print(ds)
        print(sch)