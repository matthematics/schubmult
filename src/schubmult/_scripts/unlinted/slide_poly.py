from schubmult import *
from schubmult.symbolic import S, expand_seq


def monomial_slide(comp, genset):
    if len(comp) == 0 or all(c == 0 for c in comp):
        return S.One
    
    if all(c > 0 for c in comp):
        return expand_seq(comp, genset)
    
    ret = S.Zero

    min_nonzero = min(i for i, c in enumerate(comp) if c > 0)

    power = comp[min_nonzero]
    

    for i in range(min_nonzero + 1):
        ret += genset[min_nonzero + 1 - i]**power * monomial_slide((0,) * i + comp[min_nonzero+1:], genset[min_nonzero + 1 - i:])
    
    return ret

def all_weak_comps(length, max_degree):
    if length == 0:
        yield ()
    elif length == 1:
        for i in range(max_degree + 1):
            yield (i,)
    else:
        for i in range(max_degree + 1):
            for tail in all_weak_comps(length - 1, max_degree - i):
                yield (i,) + tail

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    deg = int(sys.argv[2])
    for comp in all_weak_comps(n, deg):
        print(f"Comp: {comp}, Monomial slide: {monomial_slide(comp, Sx.genset)}")