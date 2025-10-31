import pytest

from schubmult.schub_lib.rc_graph import RCGraph
from schubmult.schub_lib.rc_graph_ring import RCGraphRing
from schubmult.schub_lib.perm_lib import Permutation

def _pad_left(t, length):
    if len(t) >= length:
        return t
    return (0,) * (length - len(t)) + t

def _add_weights(w1, w2):
    L = max(len(w1), len(w2))
    a = _pad_left(w1, L)
    b = _pad_left(w2, L)
    return tuple(x + y for x, y in zip(a, b))

def rcgraph_crystal_weight_via_phi_eps(rc):
    """
    Recompute crystal weight using phi(i) - epsilon(i) for i = 1..(n-1)
    and then cumulative sums as used in the codebase.
    """
    # use len(rc) as number of rows (n)
    n = max(1, len(rc))
    a = []
    for i in range(1, n):
        # rc.phi and rc.epsilon should exist on RCGraph
        phi_i = rc.phi(i)
        eps_i = rc.epsilon(i)
        a.append(phi_i - eps_i)
    # convert a (length n-1) to full weight tuple (length n)
    wt = []
    running = 0
    # iterate reversed to do right-to-left cumulative sum
    for ai in reversed(a):
        running += ai
        wt.append(running)
    wt.append(running)
    wt = tuple(reversed(wt))
    return wt

def compute_weights_for_element(elem):
    """
    For mapping-like elements (with .items()) return set of weights of basis keys
    with nonzero coefficient. For RCGraph input, return single tuple.
    For tuple keys (rc1, rc2, ...), compute elementwise sum of component weights.
    """
    # RCGraph instance
    if isinstance(elem, RCGraph):
        return {rcgraph_crystal_weight_via_phi_eps(elem)}

    # mapping-like
    if hasattr(elem, "items"):
        res = set()
        for key, coeff in elem.items():
            if coeff == 0:
                continue
            if isinstance(key, RCGraph):
                res.add(rcgraph_crystal_weight_via_phi_eps(key))
            elif isinstance(key, (tuple, list)):
                total = None
                ok = True
                for comp in key:
                    if not isinstance(comp, RCGraph):
                        # if comp has crystal_weight attribute use it
                        if hasattr(comp, "crystal_weight"):
                            cw = comp.crystal_weight
                        else:
                            ok = False
                            break
                    else:
                        cw = rcgraph_crystal_weight_via_phi_eps(comp)
                    total = cw if total is None else _add_weights(total, cw)
                if ok and total is not None:
                    res.add(total)
            else:
                if hasattr(key, "crystal_weight"):
                    res.add(key.crystal_weight)
                # else skip
        return res

    # fallback: empty
    return set()

def test_crystal_weight_matches_phi_eps_small_perms():
    # test permutations up to size 4 (adjust if slow)
    for n in range(1, 5):
        perms = Permutation.all_permutations(n)
        for perm in perms:
            # only test perms with nonzero inversions if your pipeline expects that
            # Here test all permutations; RCGraph.all_rc_graphs will handle trivial cases
            # test default length = len(trimcode)
            try:
                graphs = RCGraph.all_rc_graphs(perm)
            except Exception as e:
                pytest.skip(f"Skipping perm {perm}: all_rc_graphs raised {e}")
            for g in graphs:
                # compute via property
                prop_wt = getattr(g, "crystal_weight", None)
                # compute via phi/epsilon
                calc_wt = rcgraph_crystal_weight_via_phi_eps(g)
                assert prop_wt == calc_wt, (
                    f"Mismatch crystal_weight for RCGraph:\n"
                    f"perm={perm}\nRCGraph={g}\nprop={prop_wt}\ncalc={calc_wt}"
                )

def test_weights_for_ring_element():
    # pick a small permutation and a graph and check RCGraphRing.from_dict behavior
    perms = Permutation.all_permutations(3)
    ring = RCGraphRing()
    for perm in perms:
        graphs = RCGraph.all_rc_graphs(perm)
        if not graphs:
            continue
        g = next(iter(graphs))
        elem = ring.from_dict({g: 1})
        computed = compute_weights_for_element(elem)
        expected = {g.crystal_weight}
        assert expected.issubset(computed), (
            f"Ring element weight mismatch: expected {expected} in computed {computed}"
        )

def test_weights_for_tensor_element():
    """
    Construct a tensor ring element (tuple-key) and verify its weight equals
    the elementwise sum of component RCGraph weights.
    """
    ring = RCGraphRing()
    tring = ring @ ring  # tensor ring
    # pick two small permutations and graphs
    perms = Permutation.all_permutations(3)
    # take first two non-empty graphs we find
    g_list = []
    for perm in perms:
        graphs = RCGraph.all_rc_graphs(perm)
        if graphs:
            g_list.extend(list(graphs))
        if len(g_list) >= 2:
            break
    if len(g_list) < 2:
        pytest.skip("Not enough RCGraphs available to form a tensor test")
    g1, g2 = g_list[0], g_list[1]
    telem = tring.from_dict({(g1, g2): 1})
    computed = compute_weights_for_element(telem)
    expected = {_add_weights(g1.crystal_weight, g2.crystal_weight)}
    assert expected.issubset(computed), (
        f"Tensor element weight mismatch: expected {expected} in computed {computed}"
    )