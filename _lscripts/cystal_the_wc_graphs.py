from schubmult import *
import numpy as np
from schubmult.combinatorics.pipe_dream import PipeDream
from schubmult.symbolic.common_polys import schub_elem_sym_to_groth_elem_sym_dict, schub_dict_to_groth_dict
from schubmult.symbolic import *
from functools import cache
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
    length = n + 2
    perms = Permutation.all_permutations(n)
    interestin_perms = [perm for perm in perms if perm.inv != 0]
    for perm in perms:
        poly = untagged_groth_elem(perm, length)
        mapping = {}
        unmapping = {}
        # correct_weights = (perm * Permutation.w0(n)).pad_code(n - 1)
        weight = tuple(reversed([n - 1 - j - w for j, w in enumerate((perm * Permutation.w0(n)).pad_code(n - 1))]))
        #cw_dict = {n + 1 - i: ww for i, ww in enumerate(weight, start=1)}
        for key, coeff in poly.items():
            wc = bw.key_to_wc_graph(key).resize(n - 1)
            #wcc = bw.key_to_wc_graph(bw.make_key(k, key.size)).resize(n - 1)
            unmapping[CrystalGraphTensor(*key)] = wc
            if wc not in mapping:
                mapping[wc] = key
                
        
        hw_set = set()
        for wc in WCGraph.all_wc_graphs(perm, n - 1):
            if mapping[wc].is_highest_weight:
                hw_set.add(wc)

        crystals = {}
        for wc_hw in hw_set:
            crystals[wc_hw] = set()
            for paint in mapping[wc_hw].full_crystal:
                crystals[wc_hw].add(unmapping[CrystalGraphTensor(*paint)])

        pantalones = 0
        stt = set()
        r = GrothendieckRing(genset=Gx.genset, beta=1)
        for hw, crystal in crystals.items():
            stt.update(crystal)
            result = sum([wc.polyvalue(Gx.genset, beta=1) for wc in crystal])
            pato = KeyPoly.from_expr(result, length=n-1)
            assert len(tuple(key for key, v in pato.items() if sympify(v).expand() != 0)) == 1, f"More than one key label for crystal {crystal}"
            pantalones += pato

        assert r.from_expr(pantalones.expand()).almosteq(r(perm)), f"Grothendieck polynomial for {perm} does not match sum of crystals, r.from_expr({pantalones.expand()}) = {r.from_expr(pantalones.expand())} != r({perm}) = {r(perm)}"
        assert stt == WCGraph.all_wc_graphs(perm, n - 1), f"Crystals do not cover all WC graphs for {perm}, missing {set(WCGraph.all_wc_graphs(perm, n - 1)) - stt}"

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    main(n)