from schubmult import *
import numpy as np
from schubmult.combinatorics.pipe_dream import PipeDream
from schubmult.symbolic.common_polys import schub_elem_sym_to_groth_elem_sym_dict, schub_dict_to_groth_dict
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

def _to_wc(key):
    return bw.key_to_wc_graph(key)

def main(n):
    from schubmult.rings.polynomial_algebra import KeyPolyBasis, PolynomialAlgebra
    KeyPoly = PolynomialAlgebra(KeyPolyBasis(Gx.genset))
    length = n + 2
    Gx1 = GrothendieckRing(Gx.genset, beta=1)
    for perm in Permutation.all_permutations(n):
        print("Myperm is ", perm)
        tagged0 = untagged_groth_elem(perm, length)
        if perm.inv == 0:
            continue
        tagged = tagged0.to_wc_graph_ring_element()
        # assert all((v.subs(Gx._beta, 1) == 1 and wc.perm == perm) for wc, v in tagged.items() if v.subs(Gx._beta, 1) != 0), f"WCGraph mismatch for {perm.trimcode} in S_{n}, {perm=}, {[(wc.perm, v) for wc, v in tagged.items() if v != 0]}"
        
        crystals = set()
        #hw_stinkbat = {}
        for key, coeff in tagged0.items():
            # if bw.key_to_wc_graph(key).perm == perm:
            #     crys = frozenset([bw.key_to_wc_graph(kk) for kk in key.full_squared_crystal if bw.key_to_wc_graph(kk).perm == perm])
            #     crystals[crys] = 1#coeff.subs(Gx._beta, 1)
            #squared_key = CrystalGraphTensor(*(kk.square for kk in key.factors))#.to_highest_weight()[0]
            #hw_key = bw.make_key([kk.unwrap() for kk in squared_key.factors], key.size)
            #if squared_key.is_highest_weight and bw.key_to_wc_graph(key).perm == perm:
            # if key.is_highest_weight and bw.key_to_wc_graph(key).perm == perm:
            #     hw_stinkbat[key] = 1#coeff.subs(Gx._beta, 1)
            if any(key in c for c in crystals):
                continue
            if _to_wc(key).perm == perm:
                the_set = frozenset(key.full_crystal_bothways(lambda x: _to_wc(x).perm == perm and _to_wc(x).hecke_invariant[0] == _to_wc(key).hecke_invariant[0]))
                if the_set:
                    crystals.add(the_set)
        
        poly = 0
        #for crys, coeff in hw_stinkbat.items():
        for crys in crystals:
            spanko = 0
            # squared_key = CrystalGraphTensor(*(kk.square for kk in crys.factors))
            #wc_base = bw.key_to_wc_graph(crys)
            # for keykeykey in squared_key.full_crystal_bothways:
            #     wc = bw.key_to_wc_graph(bw.make_key([kk.unwrap() for kk in keykeykey.factors], crys.size))
            #     if wc.perm != perm:
            #         continue
            # for keykeykey in crys.full_crystal_bothways:
            #     wc = bw.key_to_wc_graph(keykeykey)
            #     if wc.perm != perm:
            #         continue
            #     spanko += coeff * wc.polyvalue(Gx.genset, beta=1)
            for key in crys:
                #wc = bw.key_to_wc_graph(key)
                spanko += _to_wc(key).polyvalue(Gx.genset, beta=1)
                #spanko += coeff *    keykey[0].polyvalue(Gx.genset, beta=1)
            print(KeyPoly.from_expr(spanko))
        
            poly += spanko
        poly = poly.expand()
        print("Num crystals for ", perm, " is ", len(crystals))
        # if len(crystals) > 1:
        #     input()
        assert (poly - Gx1(perm).expand()).expand() == 0, f"Failed for {perm.trimcode} in S_{n}, got {poly} vs {Gx1(perm).expand()}"
        

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    main(n)