from schubmult import *
from schubmult.symbolic import *
from schubmult.abc import *
from multiprocessing import Pool, cpu_count
from functools import partial

def is_decomposable(w):
    # return False
    for i in range(1, len(w) - 1):
        coset, w_J = w.coset_decomp(*list(range(1, i + 1)), *list(range(i + 2, len(w))))
        if coset.inv == 0 and set(w_J.code[: i + 1]) != {0} and set(w_J.code[i + 2 :]) != {0}:
            return True
    return False

def decompose(w):
    for i in range(1, len(w) - 1):
        coset, w_J = w.coset_decomp(*list(range(1, i + 1)), *list(range(i + 2, len(w))))
        if coset.inv == 0 and set(w_J.code[: i + 1]) != {0} and set(w_J.code[i + 2 :]) != {0}:
            return (uncode(w_J.code[: i + 1]), *decompose(uncode(w_J.code[i + 2 :])))
    return (w,)

def test_perm_pair(dom_perm, perm):
    """Test a single (dominant_perm, perm) pair."""
    
    test_prod = DSx(dom_perm) * DSx(perm, "z")
    failures = []
    
    for w, v in test_prod.items():
        result = S.Zero
        # if is_decomposable(perm):
        #     factors = decompose(perm)
        #     sponge = DSx(dom_perm)
        #     for factor in factors:
        #         sponge *= Sx(factor)
        #     for rc, coeff in sponge.items():
        #         result += coeff * Sx(rc.polyvalue(y[:len(rc)])) * expand_seq(rc.length_vector, y[len(rc):])
        #     sputnik = expand(v - result)
        #     if sputnik != S.Zero:
        #         failures.append((dom_perm, perm, w, sputnik, result, v))
        #     continue
        for rc in RCGraph.all_rc_graphs(perm, len(perm.trimcode)):
            dpset = rc.dualpieri(dom_perm, w)
            for dp in dpset:
                vlist = dp[0]
                the_rc = dp[1]
                result += prod([prod([y[i+1] - z[a] for a in vc]) for i, vc in enumerate(vlist)]) * the_rc.polyvalue(y[len(vlist):], z)

        sputnik = expand(v - result)
        if sputnik != S.Zero:
            failures.append((dom_perm, perm, w, sputnik, result, v))
    
    return failures

def process_dom_perm(dom_perm, perms):
    """Process all perms for a given dominant permutation."""
    all_failures = []
    for perm in perms:
        failures = test_perm_pair(dom_perm, perm)
        all_failures.extend(failures)
    return all_failures

if __name__ == "__main__":
    import sys

    n = int(sys.argv[1])
    num_processes = int(sys.argv[2]) if len(sys.argv) > 2 else cpu_count()
    
    # print("STNAKBAT NONDECOMP")
    perms = [perm for perm in Permutation.all_permutations(n)]
    dom_perms = [p for p in perms if p == p.minimal_dominant_above()]
    
    print(f"Testing {len(dom_perms)} dominant permutations against {len(perms)} permutations using {num_processes} processes")
    
    # Parallelize over dominant permutations
    with Pool(num_processes) as pool:
        worker = partial(process_dom_perm, perms=perms)
        results = pool.map(worker, dom_perms)
    
    # Collect all failures
    all_failures = []
    for failures in results:
        all_failures.extend(failures)
    
    if all_failures:
        print(f"\n{len(all_failures)} test(s) failed:")
        for dom_perm, perm, w, sputnik, result, v in all_failures:
            print(f"Failed for {dom_perm} and {perm} with difference {sputnik} {result=} {w=} {v=}")
        sys.exit(1)
    else:
        print(f"\nAll tests passed!")
        sys.exit(0)