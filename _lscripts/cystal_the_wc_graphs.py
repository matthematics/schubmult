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
def untagged_groth_elem(perm, length, betasub=1):
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

    result = bw.from_dict({k: v.subs(Gx._beta, betasub) for k, v in result.items()})
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
            # if wc not in mapping:
            spinach = key.raising_operator
            spinach2 = key.lowering_operator
            wc.raising_operator = lambda i: bw.key_to_wc_graph(bw.make_key(spinach(i), key.size)) if spinach(i) is not None else None
            wc.lowering_operator = lambda i: bw.key_to_wc_graph(bw.make_key(spinach2(i), key.size)) if spinach2(i) is not None else None

            mapping[wc] = key
            
        
        hw_set = set()
        for wc in WCGraph.all_wc_graphs(perm, n - 1):
            if wc.is_highest_weight:
                hw_set.add(wc)

        crystals = {}
        for wc_hw in hw_set:
            crystals[wc_hw] = set()
            for wcc in wc_hw.full_crystal:
                crystals[wc_hw].add(wcc)

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

def test_crystal_well_defined(n, verbose=True):
    """Test whether the induced crystal structure on WCGraphs is well-defined.

    For each permutation, ``untagged_groth_elem(perm, length)`` produces a set of
    tensor keys.  Each key maps to a WCGraph via ``key_to_wc_graph(key).resize(n-1)``,
    but this map is many-to-one: several keys can represent the same WCGraph.

    We define the crystal operators on a WCGraph by picking a representative key and
    applying the tensor operator, then mapping back:

        wc.raising_operator(i)  := key_to_wc_graph(rep_key.raising_operator(i))
        wc.lowering_operator(i) := key_to_wc_graph(rep_key.lowering_operator(i))

    This is only well-defined if *every* representative of a given WCGraph yields the
    same image under each operator.  This routine checks exactly that and reports any
    WCGraph/index/operator where representatives disagree.
    """
    length = n + 2
    max_index = length
    perms = Permutation.all_permutations(n)

    def apply(key, op_name, i):
        op = getattr(key, op_name)
        res = op(i)
        if res is None:
            return None
        return bw.key_to_wc_graph(bw.make_key(res, key.size)).resize(n - 1)

    total_classes = 0
    total_multi = 0
    failures = []
    for perm in perms:
        if perm.inv == 0:
            continue
        poly = untagged_groth_elem(perm, length)

        # Group representative keys by the WCGraph they map to.
        groups = {}
        for key in poly:
            wc = bw.key_to_wc_graph(key).resize(n - 1)
            groups.setdefault(wc, []).append(key)

        total_classes += len(groups)
        for wc, keys in groups.items():
            if len(keys) < 2:
                continue
            total_multi += 1
            for i in range(1, max_index + 1):
                for op_name in ("raising_operator", "lowering_operator"):
                    images = set()
                    has_none = False
                    for key in keys:
                        img = apply(key, op_name, i)
                        if img is None:
                            has_none = True
                        else:
                            images.add(img)
                    distinct = len(images) + (1 if has_none else 0)
                    if distinct > 1:
                        failures.append((perm, wc, op_name, i, len(keys), len(images), has_none))
                        if verbose:
                            print(
                                f"[FAIL] perm={perm} op={op_name} i={i}: "
                                f"{len(keys)} reps -> {len(images)} distinct WCGraphs"
                                f"{' + None' if has_none else ''}\n       wc={wc}",
                            )

    print(
        f"n={n}: {total_classes} WCGraph classes, {total_multi} with multiple key reps, "
        f"{len(failures)} well-definedness failures.",
    )
    if not failures:
        print("Induced crystal structure is WELL-DEFINED across all representatives.")
    else:
        print("Induced crystal structure is NOT well-defined (representatives disagree).")
    return failures


def test_wc_to_rc_rep_independent(n, verbose=True):
    """Check that the WCGraph -> PRINCIPAL RCGraph map is well-defined.

    ``untagged_groth_elem(perm, length)`` produces many tensor keys; several distinct
    keys can evaluate to the SAME WCGraph via ``bw.key_to_wc_graph(key)``.  Grouping
    by the WCGraph alone (or ``(wc, degree_tuple)``, or restricting to highest-weight
    keys) is NOT well-defined (all fail at n=5).

    Here we convert each key factor-wise (``_convert_elem_rc`` then
    ``br.key_to_rc_graph``) and keep ONLY the resulting RCGraphs that are PRINCIPAL
    (``rc.is_principal``, i.e. ``rc.perm == uncode(rc.length_vector)``).  We then ask
    whether each evaluated WCGraph yields a UNIQUE principal RCGraph.
    """
    length = n + 2
    perms = Permutation.all_permutations(n)

    def key_to_rc(key):
        rc_factors = tuple(wc._convert_elem_rc() for wc in key)
        return br.key_to_rc_graph(br.make_key(rc_factors, key.size))

    total_classes = 0
    total_multi = 0
    failures = []
    crashes = []
    for perm in perms:
        if perm.inv == 0:
            continue
        poly = untagged_groth_elem(perm, length)

        # Convert every key; group principal RCGraph images by evaluated WCGraph.
        groups = {}
        for key in poly:
            wc = bw.key_to_wc_graph(key).resize(n - 1)
            try:
                rc = key_to_rc(key)
            except Exception as exc:  # noqa: BLE001
                crashes.append((perm, wc, key, repr(exc)))
                print(f"[CRASH] perm={perm}: key_to_rc failed")
                print(f"        wc={wc}")
                print(f"        key factors (deg,desc,wt): {[(len(w.perm_word), w.perm.max_descent, w.length_vector) for w in key]}")
                print(f"        error: {exc}")
                continue
            if not rc.is_principal:
                continue
            groups.setdefault(wc, set()).add(rc)

        total_classes += len(groups)
        for wc, images in groups.items():
            if len(images) < 2:
                if verbose:
                    (rc,) = tuple(images)
                    print(f"[ok] perm={perm}: unique principal RCGraph (perm {rc.perm})")
                continue
            total_multi += 1
            failures.append((perm, wc, len(images)))
            print(f"[FAIL] perm={perm}: WCGraph -> {len(images)} distinct PRINCIPAL RCGraphs")
            print(f"       wc={wc}")
            for rc in images:
                print(f"       -> rc.perm={rc.perm}\n{rc}")
            raise AssertionError(f"WCGraph->principal RCGraph map not unique for perm={perm}, wc={wc}")

    print(
        f"n={n}: {total_classes} WCGraph classes with a principal RCGraph image, "
        f"{len(failures)} non-unique, {len(crashes)} conversion crashes.",
    )
    if not failures and not crashes:
        print("WCGraph -> RCGraph map is UNIQUE when restricted to principal RCGraph images.")
    elif failures:
        print("WCGraph -> principal RCGraph map is NOT unique.")
    else:
        print("No uniqueness failures, but some keys could not be converted.")
    return failures


def test_wc_by_degree_and_rcperm(n, verbose=True):
    """Experimental: injective key -> WCGraph among canonically-factored keys.

    We start with the keys from the Grothendieck representation
    ``untagged_groth_elem(perm, length)``.  For a key ``K`` in ``bw``:

    - ``wc = bw.key_to_wc_graph(K)``             (the WCGraph)
    - ``rc_key = convert(K)``                    (br key of RCGraph elem syms; same
      degrees/descents as ``K`` by construction)
    - ``rc = br.key_to_rc_graph(rc_key)``        (the end-result RCGraph)

    FINAL condition (the one that fixes the degree/descent signature, i.e. the
    permutation has coefficient 1 in the complete product): keep ``K`` ONLY when

        rc_key == br.from_rc_graph(rc, size)

    i.e. ``rc_key`` is exactly the canonical factorization of its own end-result
    RCGraph.  Among the keys passing this filter, group by ``(signature, wc)`` where
    ``signature = tuple((len(f.perm_word), f.perm.max_descent) for f in K)``, and check
    that each ``(signature, wc)`` corresponds to a UNIQUE key -- i.e. no other
    canonically-factored key with the same degrees/descents maps to the same WCGraph.
    """
    length = n + 2
    perms = Permutation.all_permutations(n)

    def convert(key):
        rc_factors = tuple(wc._convert_elem_rc() for wc in key)
        return br.make_key(rc_factors, key.size)

    def signature(key):
        return tuple((len(f.perm_word), f.perm.max_descent) for f in key)

    # (signature, wc) -> set of canonically-factored bw keys, global.
    index = {}
    canonical = 0
    non_canonical = 0
    crashes = 0
    for perm in perms:
        if perm.inv == 0:
            continue
        poly = untagged_groth_elem(perm, length)
        for key in poly:
            try:
                wc = bw.key_to_wc_graph(key)
                rc_key = convert(key)
                rc = br.key_to_rc_graph(rc_key)
                (recovered_rc_key,) = tuple(br.from_rc_graph(rc, key.size))
            except Exception:  # noqa: BLE001
                crashes += 1
                continue
            # FINAL condition: rc_key must be the canonical factorization of its RCGraph.
            if recovered_rc_key != rc_key:
                non_canonical += 1
                continue
            canonical += 1
            index.setdefault((signature(key), wc), set()).add(key)

    failures = []
    for (sig, wc), keys in index.items():
        if len(keys) > 1:
            failures.append((sig, wc, len(keys)))
            print(f"[FAIL] signature={sig}: {len(keys)} canonical keys with this (degree,descent) map to the same WCGraph")
            print(f"       wc={wc}")
            for key in keys:
                print(f"       key={key}")
            raise AssertionError(f"key -> WCGraph not injective within signature {sig} (wc={wc})")

    print(
        f"n={n}: {canonical} canonically-factored keys ({non_canonical} non-canonical, {crashes} crashes), "
        f"{len(index)} (signature, WCGraph) classes, {len(failures)} with multiple keys.",
    )
    if not failures:
        print("Among canonically-factored keys, key -> WCGraph is INJECTIVE within a fixed (degree, descent) signature.")
    else:
        print("key -> WCGraph is NOT injective even among canonically-factored keys of a fixed signature.")
    return failures


def test_groth_to_schub_via_rc(n, verbose=True):
    """Change basis Grothendieck -> Schubert through the untagged tensor representation.

    For each permutation:

    1. ``untagged_groth_elem(perm, length)`` gives a ``bw`` element whose keys are
       tensors of WCGraph (Grothendieck) elem syms, coefficients already at beta=1.
    2. Each WCGraph elem-sym factor is replaced by the RCGraph elem sym of the SAME
       weight and SAME descent (``RCGraph.elem_sym_rcs(len(perm_word), perm.max_descent,
       weight=length_vector)``, i.e. ``WCGraph._convert_elem_rc()``); the degree bumps
       from the reduced degree up to ``len(perm_word)``.  This turns each ``bw`` key
       into a ``br`` (BoundedRCFactorAlgebra) key of the same size.
    3. ``to_rc_graph_ring_element()`` evaluates the ``br`` element to a sum of RCGraphs.
    4. Collect ``{rc.perm: coefficient}`` and compare against
       ``WCGraph.groth_to_schub(perm, beta=1)``.
    """
    length = n + 2
    perms = Permutation.all_permutations(n)

    failures = []
    tested = 0
    for perm in perms:
        if perm.inv == 0:
            continue
        tested += 1

        poly = untagged_groth_elem(perm, length)

        # Rebuild each bw key as a br key of RCGraph elem syms (same weight, same descent).
        br_elem = br.zero
        for key, coeff in poly.items():
            rc_factors = tuple(wc._convert_elem_rc() for wc in key)
            br_key = br.make_key(rc_factors, key.size)
            br_elem += coeff * br(br_key)

        # Evaluate to RCGraphs and collect perm -> coefficient.
        # Different RCGraphs are parallel copies: do NOT sum their coefficients.
        # RCGraphs sharing a perm should carry the same coefficient; take it (and
        # flag any disagreement rather than accumulating).
        rc_ring_elem = br_elem.to_rc_graph_ring_element()
        got = {}
        parallel_conflict = {}
        for rc, coeff in rc_ring_elem.items():
            c = sympify(coeff).expand()
            if c == 0:
                continue
            if rc.perm in got and got[rc.perm] != c:
                parallel_conflict[rc.perm] = (got[rc.perm], c)
            got[rc.perm] = c
        got = {p: c for p, c in got.items() if sympify(c).expand() != 0}

        # Reference change of basis at beta=1.
        want = {p: sympify(c).subs(Gx._beta, 1) for p, c in WCGraph.groth_to_schub(perm, Gx._beta).items()}
        want = {p: sympify(c).expand() for p, c in want.items() if sympify(c).expand() != 0}

        if got != want:
            failures.append(perm)
            mismatch_got = {p: got[p] for p in set(got) | set(want) if got.get(p, 0) != want.get(p, 0)}
            mismatch_want = {p: want.get(p, 0) for p in set(got) | set(want) if got.get(p, 0) != want.get(p, 0)}
            print(f"[FAIL] perm={perm}")
            print(f"       via RC:         {got}")
            print(f"       groth_to_schub: {want}")
            print(f"       mismatch got:  {mismatch_got}")
            print(f"       mismatch want: {mismatch_want}")
            if parallel_conflict:
                print(f"       PARALLEL CONFLICT (same perm, different coeff): {parallel_conflict}")
            raise AssertionError(f"Groth->Schub via RC mismatch for perm={perm}: via RC {got} != groth_to_schub {want}")
        if verbose:
            print(f"[ok] perm={perm}: {len(got)} Schubert terms match")

    print(f"n={n}: {tested} permutations tested, {len(failures)} mismatches.")
    if not failures:
        print("Groth->Schub via RCGraph elem-sym representation AGREES with groth_to_schub (beta=1).")
    else:
        print(f"Mismatches for: {failures}")
    return failures


if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    mode = sys.argv[2] if len(sys.argv) > 2 else "main"
    if mode == "crystal":
        test_crystal_well_defined(n)
    elif mode == "g2s":
        test_groth_to_schub_via_rc(n)
    elif mode == "wc2rc":
        test_wc_to_rc_rep_independent(n)
    elif mode == "wcdeg":
        test_wc_by_degree_and_rcperm(n)
    else:
        main(n)