from schubmult import *
import numpy as np
from schubmult.combinatorics.indexed_forests import grove_polynomial
from schubmult.combinatorics.pipe_dream import PipeDream
from schubmult.symbolic.common_polys import schub_elem_sym_to_groth_elem_sym_dict, schub_dict_to_groth_dict
from schubmult.symbolic import *
from functools import cache
br = BoundedRCFactorAlgebra()
bw = BoundedWCFactorAlgebra()
from sympy import pretty_print, sympify, Pow, Mul
from schubmult.rings.polynomial_algebra import Grove, GrovePolyBasis

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


@cache
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

    result = bw.from_dict({k: v.subs(Gx._beta, 1) for k, v in result.items()})
    result = bw.from_dict({k: v for k, v in result.items() if v != 0})
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
    from schubmult.rings.polynomial_algebra import KeyPolyBasis, PolynomialAlgebra
    import itertools
    KeyPoly = PolynomialAlgebra(KeyPolyBasis(Gx.genset))
    length = n - 1
    perms = Permutation.all_permutations(n)
    interestin_perms = [perm for perm in perms if perm.inv != 0]
    wcr = WCGraphRing()
    Grove1 = PolynomialAlgebra(GrovePolyBasis(Gx.genset, beta=1))

    # --- Memoization of per-perm / per-key work ---------------------------
    # Each perm appears in O(#perms) combinations, so its per-perm work is
    # otherwise recomputed many times. These caches keep the exact same
    # computation path, just evaluated once per distinct argument.
    @cache
    def _key_to_wc(k):
        # squash_product bijection; expensive and reused across many perms
        return bw.key_to_wc_graph(k)

    @cache
    def _grove_expand(comp):
        return Grove1(*comp).expand()

    @cache
    def _gx_of(comp):
        # Gx.from_expr is a ring homomorphism, so caching this per single comp
        # and multiplying in the Gx ring reproduces Gx.from_expr(prod of exprs)
        # without the per-pair symbolic polynomial multiply + re-parse.
        return Gx.from_expr(_grove_expand(comp))

    @cache
    def _all_wcs(perm):
        return WCGraph.all_wc_graphs(perm, n - 1)

    @cache
    def _grove_elem(weight):
        return Grove1(*weight)

    @cache
    def _grove_weights_of(k):
        # Forest weights of the WCGraphs of k that survive the grove filter.
        # In the real_product loop every such wc carries the same coefficient,
        # so we only need the (multiset of) qualifying weights.
        return tuple(wc.forest_weight for wc in _all_wcs(k) if wc.forest_weight == wc.length_vector)

    @cache
    def _poly_for(perm):
        comp = perm.pad_code(n - 1)
        gwc = WCGraph.grove_wcs(comp)

        def _keep(k):
            wcg = _key_to_wc(k)
            return wcg.resize(n - 1) in gwc or wcg.perm != perm

        return bw.from_dict({k: v for k, v in untagged_groth_elem(perm, length).items() if _keep(k)})

    for perm1, perm2 in itertools.combinations(interestin_perms, 2):
        # perm1 = uncode([0,0,0,1])
        # perm2 = uncode([0,0,2,1])
        comp1 = perm1.pad_code(n-1)
        comp2 = perm2.pad_code(n-1)
        poly1 = _poly_for(perm1)
        poly2 = _poly_for(perm2)
        if any(v < 0 for v in poly1.values()) or any(v < 0 for v in poly2.values()):
            print(f"Negative coefficient")# for {comp1=}, {comp2=}: {poly1=}, {poly2=}")
            
        try_poly = (poly1 * poly2).to_wc_graph_ring_element()
        
        real_product_g = Gx.from_dict({k: sympify(v).subs(Gx._beta, 1) for k, v in (_gx_of(comp1) * _gx_of(comp2)).items()})
        real_product = 0
        for k, v in real_product_g.items():
            for fw in _grove_weights_of(k):
                real_product += v * _grove_elem(fw)
            
        painted_potato = 0
        for wc, v in try_poly.items():
            #assert sympify(real_product.get(wc.perm, 0)).subs(Gx._beta, 1) == v, f"Mismatch for {perm1=}, {perm2=}: {try_poly=} vs {real_product=}"
            if wc.forest_weight == wc.length_vector:
                painted_potato += v * _grove_elem(wc.forest_weight)
        assert real_product.almosteq(painted_potato), f"Mismatch for {perm1=}, {perm2=}: {(painted_potato-real_product)=}"
        print(f"Success {comp1=}, {comp2=}")
        # break
        # tagged = tagged0.to_wc_graph_ring_element()
        # assert all((v == 1 and wc.perm == perm) for wc, v in tagged.items()), f"WCGraph mismatch for {perm.trimcode} in S_{n}, {perm=}, {[(wc.perm, v) for wc, v in tagged.items() if v != 0]}"
        
        # hw_dict = {}
        # for key, coeff in tagged.items():
        #     if key.is_highest_weight:
        #         hw_dict[key] = hw_dict.get(key, 0) + coeff.subs(Gx._beta, 1)

        # poly = 0
        # for key_hw, coeff in hw_dict.items():
        #     # properly_sortable = [rc for rc in key_hw.full_crystal if tuple(sorted(rc.crystal_weight, reverse=True)) == tuple(key_hw.crystal_weight)]
        #     # min_vec = min([(i, np.cumsum(properly_sortable[i].crystal_weight).tolist()) for i in range(len(properly_sortable))], key=lambda x: x[1])[0]
        #     spanko = 0
        #     # coeff2, mulmul = sympify(coeff).as_coeff_Mul()
        #     # assert not sympify(coeff2).is_negative, f"Negative coefficient {coeff2} for {perm.trimcode} in S_{n}"
        #     #spinach = 0
        #     for keykey in key_hw.full_crystal:
        #         spanko += coeff * keykey.polyvalue(Gx.genset, beta=Gx._beta, prop_beta=True)
        #         #spinach += coeff * keykey.polyvalue(Gx.genset, beta=1, prop_beta=False)

        #     # keyspanko = KeyPoly.from_expr(spinach)
        #     # assert len([porko for porko, v in keyspanko.items() if sympify(v).expand() != 0]) <= 1, f"Key polynomial expansion for {perm.trimcode} in S_{n} has more than one term: {keyspanko}"
        #     poly += spanko
        # poly = poly.expand()
        
        # poly = tagged.polyvalue(Gx.genset, beta=Gx._beta, prop_beta=True)
        # assert (poly - Gx(perm).expand()).expand() == 0, f"Failed for {perm.trimcode} in S_{n}, got {poly} vs {Gx(perm).expand()}"
        # #spinach = reduce_it.to_rc_graph_ring_element()
        # # reduced_potato = reduce_it.to_rc_graph_ring_element().resize(length)
        # # dct2 = WCGraph.groth_to_schub(perm, Gx._beta)
        # # for key, rc_key in dct.items():
        # #     rcc = br.key_to_rc_graph(rc_key)
        # #     if reduce_it.get(rc_key, 0) == 0:
        # #         continue
        # #     if reduced_potato.get(br.key_to_rc_graph(rc_key).resize(length), 0) == 0:
        # #         continue
        # #     excess = 0
        # #     if not(sympify(reduce_it.get(rc_key, 0)).is_Number):
        # #         coeff, popp = sympify(reduce_it.get(rc_key, 0)).as_coeff_Mul()
        # #         if popp == Gx._beta:
        # #             excess = 1
        # #         elif popp.is_Pow and popp.base == Gx._beta:
        # #             excess = popp.exp
        # #         else:
        # #             raise ValueError(f"Unexpected coefficient {reduce_it.get(rc_key, 0)} for {perm.trimcode} in S_{n}")
        # #     if len(key) == 1:
        # #         continue
        # #     # if excess + sum([kk.excess for kk in key]) != sum([len(kk.perm_word) for kk in key]) - perm.inv:
        # #     #     continue
        # #     # copi = PipeDream.from_rc_graph(rc.resize(length)).co_pipe_dream()
        # #     # test_perm = copi.perm * Permutation.w0(copi.rows)
        # #     # print("Pago")
        # #     # pretty_print(rc)
        # #     # assert test_perm == perm, f"RC graph mismatch for {perm.trimcode} in S_{n}: got {test_perm} vs {perm}"
        # #     # if br.key_to_rc_graph(rc_key).is_principal:
        # #     #     print("PRINCIPAL")
        # #     #     pretty_print(bw(bw.make_key(key, key.size)))
        # #     #assert dct.get(rc.perm, 0) == coeff, f"RC graph mismatch for {perm.trimcode} in S_{n}: got {coeff} vs {dct.get(rc.perm, 0)}\n{rc=}"
        # #     key_hw = key.to_highest_weight()[0]
        # #     lv = [0] * length
        # #     for kk in key_hw:
        # #         for i, j in enumerate(kk.length_vector):
        # #             lv[i] += j
        # #     wcs = {wcc for wcc in WCGraph.all_wc_graphs(perm, length, weight=lv) if wcc.excess <= excess + sum([kk.excess for kk in key])}
        # #     assert len(wcs) > 0, f"No WCGraphs found for {perm.trimcode} in S_{n} with weight {lv} and highest weight {key_hw}"
        # #     if len(wcs) > 1:
        # #         print("Petunia")
        # #         pretty_print(key_hw)
        # #         print("Crampy bra")
        # #         pretty_print(br.key_to_rc_graph(rc_key).normalize())
        # #         for wc in wcs:
        # #             print("Candidate")
        # #             pretty_print(wc)

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    main(n)