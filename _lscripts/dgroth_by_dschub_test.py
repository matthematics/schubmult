from schubmult import *
from schubmult.symbolic.common_polys import grothendieck_poly_with_ring
from schubmult.symbolic import S, sympify, prod, expand_func, expand, Integer, Mul, Add, Pow, Symbol
from schubmult.abc import E, y, z, H
from schubmult.symbolic.common_polys import efficient_subs
from functools import cache
from schubmult.symbolic.common_polys import _vars, efficient_subs, elem_func_func_mul, elem_sym_func, elem_sym_poly
from schubmult.symbolic.poly.variables import CustomGeneratingSet, GeneratingSet, GeneratingSet_base, MaskedGeneratingSet
from schubmult.symbolic.symmetric_polynomials import FactorialElemSym
from schubmult.utils.logging import get_logger
from schubmult.utils.perm_utils import add_perm_dict, add_perm_dict_with_coeff
from schubmult.utils.schub_lib import (
    compute_vpathdicts,
    elem_sym_perms,
    elem_sym_perms_op,
    elem_sym_positional_perms,
    pull_out_var,
)
from schubmult.mult.separated_descents import grothmult_double_plus
#from schubmult.utils._mul_utils import add_perm_dict_with_coeff
beta = Gx._beta

def _check_perm(perm_list):
    import time
    start = time.monotonic()
    perm = Permutation(perm_list)
    ring = DSx([]).ring
    _beta = Gx._beta
    groth1 = grothendieck_poly_with_ring(perm, ring=ring, beta=_beta, keep_as_schub=True)
    
    for coeff in groth1.values():
        #coeff = expand(coeff)
        if isinstance(coeff, int):
            if coeff < 0:
                return (perm_list, f"Grothendieck polynomial for {perm} has negative coefficient: {coeff=}, \n{groth1=}", time.monotonic() - start)
        else:
            result = True
            def walk_expr(args):
                nonlocal result
                for arg in args:
                    if isinstance(arg, int | Integer):
                        if arg < 0:
                            result = False
                            return
                    elif isinstance(arg, Mul | Add | Pow):
                        walk_expr(arg.args)
                    elif isinstance(arg, Symbol):
                        return
                    else:
                        raise ValueError(f"Unexpected type {type(arg)} in expression {arg}")
            walk_expr(coeff.args)
            if not result:
                return (perm_list, f"Grothendieck polynomial for {perm} has negative coefficient", time.monotonic() - start)
    elapsed = time.monotonic() - start
    return (perm_list, None, elapsed)


if __name__ == "__main__":
    import sys
    from multiprocessing import Pool, cpu_count

    n = int(sys.argv[1])
    num_procs = int(sys.argv[2]) if len(sys.argv) > 2 else cpu_count()
    perms = [perm for perm in Permutation.all_permutations(n) if perm.inv != 0]
    # hardest (highest-inversion) permutations first so a few slow stragglers
    # don't get scheduled last and stall the whole pool at the end
    perms.sort(key=lambda p: p.inv, reverse=True)
    perm_lists = [list(perm) for perm in perms]

    failures = []
    with Pool(num_procs) as pool:
        for i, (perm_list, error, elapsed) in enumerate(pool.imap_unordered(_check_perm, perm_lists), 1):
            if error is not None:
                failures.append(error)
                print(f"FAIL {perm_list}", file=sys.stderr, flush=True)
            print(f"[{i}/{len(perm_lists)}] {perm_list} {'ok' if error is None else 'FAILED'} ({elapsed:.2f}s)", flush=True)

    if failures:
        raise AssertionError("\n".join(failures))
    print(f"All {len(perm_lists)} permutations verified.", flush=True)
    