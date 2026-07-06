from schubmult import *

from schubmult.combinatorics.pipe_dream import PipeDream
from schubmult.symbolic.poly.schub_poly import schub_elem_sym_to_groth_elem_sym_dict, schub_dict_to_groth_dict
from schubmult.symbolic import *
br = BoundedRCFactorAlgebra()
bw = BoundedWCFactorAlgebra()
from sympy import pretty_print, sympify, Pow, Mul

def groth_to_schub_as_rc(groth_perm: Permutation, length):
    from schubmult.combinatorics.pipe_dream import PipeDream

    boip = WCGraph.all_wc_graphs(groth_perm, length)
    ret = 0
    for wc in boip:
        cpdb = PipeDream.from_rc_graph(wc).co_pipe_dream()
        if cpdb.is_reduced:
            w0 = Permutation.w0(cpdb.rows)
            perm = cpdb.perm * w0
            ret += (Gx._beta)**(perm.inv - groth_perm.inv) * br.schub_elem(perm, length)
    return ret



def groth_elem_to_schub_elem_as_rc(p, k, length):
    from schubmult.combinatorics.pipe_dream import PipeDream
    groth_perm = uncode([0] * (k - p) + [1] * p)
    boip = WCGraph.all_wc_graphs(groth_perm, length)
    ret = (bw @ br).zero
    for wc in boip:
        cpdb = PipeDream.from_rc_graph(wc).co_pipe_dream()
        if cpdb.is_reduced:
            w0 = Permutation.w0(cpdb.rows)
            perm = cpdb.perm * w0
            for the_rc in RCGraph.all_rc_graphs(perm, length):
                the_wc = next(iter([wc2 for wc2 in boip if wc2.perm_word == wc.perm_word and wc2.length_vector == the_rc.length_vector]))
                ret += (Gx._beta)**(perm.inv - p) * bw(bw.make_key((the_wc.resize(k),),length)) @ br(br.make_key((the_rc.resize(k),),length))
    return ret

def schub_to_groth_as_rc(schub_perm: Permutation, length):
    
    ret = (r @ rw).zero
    
    for rc in RCGraph.all_rc_graphs(schub_perm, length):
        cpdb = PipeDream.from_rc_graph(rc).co_pipe_dream()
        w0 = Permutation.w0(cpdb.rows)
        perm = cpdb.perm * w0
        for wc, _ in rw.groth(perm, length).items():
            ret += (-1)**(len(wc.perm_word) - wc.perm.inv) * r(rc) @ rw(wc)
        #ret += coeff*rw.groth(perm, length)
    return ret


def untagged_groth_elem(perm, length):
    """Express the Grothendieck polynomial of ``perm`` as a tagged tensor in ``bw @ br``.

    Pipeline (no RC graphs until the final step):

    1. Expand the Grothendieck polynomial in the Schubert basis,
       ``Sx.from_expr(Gx(perm).expand())``.
    2. Rewrite each Schubert polynomial in its CEM basis, a sum of products of
       Schubert elementary symmetrics ``e_p(x_1, ..., x_k)``.
    3. Change basis on each Schubert elem sym factor abstractly via
       :func:`schub_elem_sym_to_groth_elem_sym_dict`, so every CEM monomial becomes
       a sum of products of Grothendieck elementary symmetrics, tracked only as
       ``(p, k)`` index vectors (there is no dedicated Grothendieck-elem-sym symbol).
    4. Expand each Grothendieck-elem-sym factor ``(p, k)`` into its tagged tensor
       via :func:`groth_elem_to_schub_elem_as_rc`, whose left (``bw``) factor is the
       Grothendieck elem sym over its WCGraph monomials and whose right (``br``)
       factor is the uniquely co-pipe-dream-bijected Schubert elem sym RC graph.

    Multiplying the tagged factors in the tensor ring keeps the tag paired factor
    by factor, producing ``(gelem1 @ elem1) * (gelem2 @ elem2) * ...`` rather than a
    single top-level ``groth @ schub`` tensor.
    """
    from sympy import Add, Mul, Pow, expand, sympify

    beta = Gx._beta
    identity = bw(bw.make_key((), length))

    # Step 1: Grothendieck polynomial -> Schubert basis.
    schub_dict = Sx.from_expr(Gx(perm).expand())

    result = 0
    for schub_perm, schub_coeff in schub_dict.items():
        # Step 2: each Schubert polynomial -> CEM basis (products of Schubert elem syms).
        cem_rep = expand(sympify(Sx(schub_perm).in_CEM_basis()))

        for term in Add.make_args(cem_rep):
            scalar, rest = term.as_coeff_Mul()

            schub_factors = []
            if rest != 1:
                for factor in Mul.make_args(rest):
                    if isinstance(factor, Pow):
                        base, exponent = factor.as_base_exp()
                    else:
                        base, exponent = factor, 1
                    schub_factors.extend([base] * int(exponent))

            # Step 3: rewrite each Schubert elem sym as a sum of Grothendieck elem syms,
            # distributing the product into a sum of (coeff, [(p, k), ...]) products.
            #groth_products = [(scalar, [])]
            groth_term = identity
            for factor in schub_factors:
                conversion = schub_elem_sym_to_groth_elem_sym_dict(factor.degree, factor.numvars, beta)
                # expanded = []
                # for prod_coeff, groth_factors in groth_products:
                #     for (deg, numvars), conv_coeff in conversion.items():
                #         if deg == 0:
                #             expanded.append((prod_coeff * conv_coeff, groth_factors))
                #         else:
                #             expanded.append((prod_coeff * conv_coeff, [*groth_factors, (deg, numvars)]))
                groth_term_add = 0
                for (p, k), coeff in conversion.items():
                    # Grothendieck elem sym g_p(x_1, ..., x_k) is the column perm
                    # uncode([0]*(k - p) + [1]*p), enumerated over its k-variable WCGraphs.
                    for wc in WCGraph.all_wc_graphs(uncode([0] * (k - p) + [1] * p), k):
                        groth_term_add += coeff * bw(bw.make_key((wc,), length))
                groth_term *= groth_term_add
                #groth_products = expanded

            # Step 4: expand each Grothendieck elem sym factor into its tagged tensor.
            # for prod_coeff, groth_factors in groth_products:
            #     term_tensor = identity
            #     for deg, numvars in groth_factors:
            #         term_tensor = term_tensor * groth_elem_to_schub_elem_as_rc(deg, numvars, length)

            result += schub_coeff * scalar * groth_term

    return result



# def groth_to_schub_as_rc(groth_perm: Permutation, length):
#     ret = r.zero
#     for perm, coeff in WCGraph.groth_to_schub(groth_perm, Gx._beta).items():
#         ret += coeff*r.schub(perm, length)
#     return ret


def tensor_bensor(rw_r_elem):
    ret = (rw @ r @ rw).zero
    for (wc, rc), coeff in rw_r_elem.items():
        cpdb = PipeDream.from_rc_graph(rc).co_pipe_dream()
        w0 = Permutation.w0(cpdb.rows)
        perm = cpdb.perm * w0
        for wc2, _ in rw.groth(perm, len(rc)).items():
            ret += coeff * (-1)**(len(wc2.perm_word) - wc2.perm.inv) * rw(wc)@ r(rc) @ rw(wc2)
    return ret

def main(n):
    length = n + 2
    for perm in Permutation.all_permutations(n):
        tagged = untagged_groth_elem(perm, length)
        # print(f"Tagged Grothendieck element for {perm.trimcode} (S_{n}):")
        # pretty_print(tagged)
        # print()
        # # new_tagged = 0
        # poly = 0
        # # for (wc_fact, rc_factor), coeff in tagged.items():
        # #     new_tagged += coeff * bw(wc_fact) @ br(rc_factor).to_rc_graph_ring_element().resize(length)
        # for (wc_fact, rc), coeff in tagged.items():
        poly = sum([coeff * prod([elem.polyvalue(Sx.genset, beta=Gx._beta, prop_beta=True) for elem in key]) for key, coeff in tagged.items()])
        # groth_bucket = groth_to_schub_as_rc(perm, length)
        # poly = groth_bucket.to_rc_graph_ring_element().polyvalue(Sx.genset)
        # pretty_print(new_tagged)
        assert (poly - Gx(perm).expand()).expand() == 0, f"Failed for {perm.trimcode} in S_{n}, got {poly} vs {Gx(perm).expand()}"
        # GOOD
        for key, coeff in tagged.items():
            if len(key) == 0:
                continue
            for i in range(n - 1):
                key2 = key.raising_operator(i)
                if key2 is not None:
                    assert tagged.get(bw.make_key(key2, key.size), 0) == coeff, f"Raising operator mismatch for {perm.trimcode} in S_{n}: {key} -> {key2}, got {tagged.get(key2, 0)} vs {coeff}"
                key2 = key.lowering_operator(i) 
                if key2 is not None:
                    assert tagged.get(bw.make_key(key2, key.size), 0) == coeff, f"Lowering operator mismatch for {perm.trimcode} in S_{n}: {key} -> {key2}, got {tagged.get(key2, 0)} vs {coeff}"
        reduce_it = 0
        dct = {}
        for key, coeff in tagged.items():
            rc_key = br.make_key(tuple([key[i]._convert_elem_rc().normalize() for i in range(len(key))]), key.size)
            dct[key] = rc_key
            reduce_it += coeff * (Gx._beta**(sum([len(key[i].perm_word) for i in range(len(key))]) - sum([key[i].perm.inv for i in range(len(key))]))) * br(rc_key)
        #spinach = reduce_it.to_rc_graph_ring_element()
        reduced_potato = reduce_it.to_rc_graph_ring_element().resize(length)
        dct2 = WCGraph.groth_to_schub(perm, Gx._beta)
        for key, rc_key in dct.items():
            rcc = br.key_to_rc_graph(rc_key)
            if reduce_it.get(rc_key, 0) == 0:
                continue
            if reduced_potato.get(br.key_to_rc_graph(rc_key).resize(length), 0) == 0:
                continue
            excess = 0
            if not(sympify(reduce_it.get(rc_key, 0)).is_Number):
                coeff, popp = sympify(reduce_it.get(rc_key, 0)).as_coeff_Mul()
                if popp == Gx._beta:
                    excess = 1
                elif popp.is_Pow and popp.base == Gx._beta:
                    excess = popp.exp
                else:
                    raise ValueError(f"Unexpected coefficient {reduce_it.get(rc_key, 0)} for {perm.trimcode} in S_{n}")
            if len(key) == 1:
                continue
            # if excess + sum([kk.excess for kk in key]) != sum([len(kk.perm_word) for kk in key]) - perm.inv:
            #     continue
            # copi = PipeDream.from_rc_graph(rc.resize(length)).co_pipe_dream()
            # test_perm = copi.perm * Permutation.w0(copi.rows)
            # print("Pago")
            # pretty_print(rc)
            # assert test_perm == perm, f"RC graph mismatch for {perm.trimcode} in S_{n}: got {test_perm} vs {perm}"
            # if br.key_to_rc_graph(rc_key).is_principal:
            #     print("PRINCIPAL")
            #     pretty_print(bw(bw.make_key(key, key.size)))
            #assert dct.get(rc.perm, 0) == coeff, f"RC graph mismatch for {perm.trimcode} in S_{n}: got {coeff} vs {dct.get(rc.perm, 0)}\n{rc=}"
            key_hw = key.to_highest_weight()[0]
            lv = [0] * length
            for kk in key_hw:
                for i, j in enumerate(kk.length_vector):
                    lv[i] += j
            wcs = {wcc for wcc in WCGraph.all_wc_graphs(perm, length, weight=lv) if wcc.excess <= excess + sum([kk.excess for kk in key])}
            assert len(wcs) > 0, f"No WCGraphs found for {perm.trimcode} in S_{n} with weight {lv} and highest weight {key_hw}"
            if len(wcs) > 1:
                print("Petunia")
                pretty_print(key_hw)
                print("Crampy bra")
                pretty_print(br.key_to_rc_graph(rc_key).normalize())
                for wc in wcs:
                    print("Candidate")
                    pretty_print(wc)

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    main(n)