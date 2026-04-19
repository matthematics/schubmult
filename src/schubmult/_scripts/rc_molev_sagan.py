from schubmult import *
from schubmult.abc import *
from sympy import expand, Add, Mul, S, sympify, pretty_print

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    r = RCGraphRing()
    for perm in perms:
        print("=" * 80)
        print(f"PERMUTATION: {perm}  (inv={perm.inv}, descents={perm.descents()}, trimcode={perm.trimcode})")
        print("=" * 80)

        the_cem = expand(DSx(perm).in_SEM_basis())
        print(f"  CEM expansion (in_SEM_basis): {the_cem}")

        terms = []
        for the_term in Add.make_args(sympify(the_cem)):
            coeff, rest = the_term.as_coeff_Mul()
            if rest == 1:
                terms.append((coeff, ()))
            else:
                factors = Mul.make_args(rest)
                terms.append((coeff, sorted(factors, key=lambda x: x.numvars)))

        print(f"\n  Parsed terms ({len(terms)} total):")
        for ti, term in enumerate(terms):
            print(f"    term[{ti}]: coeff={term[0]}, num_factors={len(term[1])}")
            for fi, factor in enumerate(term[1]):
                print(f"      factor[{fi}]: {factor}")
                print(f"        type={type(factor).__name__}, degree={factor.degree}, numvars={factor.numvars}")
                print(f"        genvars={factor.genvars}")
                print(f"        coeffvars={factor.coeffvars}")

        the_result = 0
        for ti, term in enumerate(terms):
            print(f"\n  --- Processing term[{ti}]: coeff={term[0]} ---")
            base = r.one
            print(f"    base (initial): {base}")

            for fi, factor in enumerate(term[1]):
                p = factor.degree
                k = factor.numvars
                #print(f"\n    factor[{fi}]: {factor}  (p={p}, k={k})")

                perm_for_elem = uncode([0] * (k - p) + [1] * p)
                all_elem_rcs = list(RCGraph.all_rc_graphs(perm_for_elem))
                #print(f"      all_rc_graphs count: {len(all_elem_rcs)}")

                yvars = DSx([]).ring.coeff_genset[1:]
                zvars = factor.coeffvars
                # print(f"      yvars (coeff_genset[1:]): {yvars[:k+2]}...")
                # print(f"      zvars (factor.coeffvars): {zvars}")
                # print(f"      len(zvars)={len(zvars)}")

                factor_result = r.zero
                for ri, elem_sym_rc in enumerate(all_elem_rcs):
                    weight = elem_sym_rc.length_vector
                    # print(f"\n      elem_sym_rc[{ri}]: perm={elem_sym_rc.perm}, len={len(elem_sym_rc)}")
                    # print(f"        weight={weight}")
                    # print(f"        rows: {[tuple(row) for row in elem_sym_rc]}")

                    try:
                        contribution = base.resize(k).double_elem_sym_squash(weight, yvars, zvars)
                        # print(f"        contribution keys ({len(list(contribution.keys()))} total): {list(contribution.keys())}")
                        # print(f"        contribution: {contribution}")
                        factor_result += contribution
                    except Exception as exc:
                        print(f"        *** EXCEPTION in double_elem_sym_squash: {type(exc).__name__}: {exc}")
                        import traceback
                        traceback.print_exc()
                        raise

                base = factor_result
                # print(f"      factor_result keys ({len(list(base.keys()))} total): {list(base.keys())}")
                # print(f"      factor_result: {base}")

            the_result += term[0] * base
            # print(f"    term[{ti}] contribution: {term[0]} * {base}")

        print(f"\n  FINAL RESULT for {perm}:")
        pretty_print(the_result)
        print(f"\n  CEM expansion was:")
        pretty_print(the_cem)
        print()