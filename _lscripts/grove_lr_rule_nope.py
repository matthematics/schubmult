from schubmult import *
from schubmult.combinatorics.indexed_forests import *
from schubmult.rings.combinatorial.forest_rc_ring import ForestRCGraphRing
from schubmult.rings.polynomial_algebra import *
from schubmult.rings.free_algebra import *
from schubmult.symbolic.common_polys import *
from schubmult.symbolic.poly.variables import *
from schubmult.utils._mul_utils import add_perm_dict
from schubmult.utils.tuple_utils import pad_tuple

import itertools

f = ForestRCGraphRing()

T = ThompsonAlgebra()
br = BoundedRCFactorAlgebra()

def _g_operator(rc, index, beta):
    """The operator G_i on rc graphs, which is the composition of the quasi-shift and trim-descent operations."""
    return f(rc).forest_trim(index) - beta * f(rc).quasi_shift(index)

def _g_grove_extractor(rc_elem, indexes, beta):
    """The operator G_i on rc graphs, which is the composition of the quasi-shift and trim-descent operations."""
    result = f.zero
    if len(indexes) > 0:
        desc, is_left_child = indexes[-1]      
        for rc, coeff in rc_elem.items():
            forest = weak_composition_to_indfor(rc.forest_weight)
            if desc not in forest.trim_descents:
                continue
            if is_left_child:
                result += coeff * (_g_operator(rc, desc, beta) + beta * f(rc).quasi_shift(desc))
            else:
                result += coeff * (_g_operator(rc, desc, beta) + beta * f(rc).quasi_shift(desc + 1))
        if len(indexes) > 1:
            return _g_grove_extractor(result, indexes=indexes[:-1], beta=beta)
    else:
        result = rc_elem
    result_values = [v for k, v in result.items() if k.perm.inv == 0]
    return sum(result_values)

# def grove_as_thompson(comp, beta, length):
#     forest = weak_composition_to_indfor(comp)
#     desc = forest.trim_descents[0]
#     is_left_child = forest.is_left_child(desc)
#     trimmed_forest = forest.trim_descent(desc)
#     return T((desc,)) * ( grove_as_thompson(trimmed_forest.code, beta, length) 
def act_on_grove_dict(dct, beta, length):
    result = {}
    for comp, coeff in dct.items():
        forest = weak_composition_to_indfor(comp)
        if len(forest.trim_descents) == 0:
            result[comp] = result.get(comp, 0) + coeff
            continue
        desc = forest.trim_descents[0]
        trimmed_forest = forest.trim_descent(desc)
        trimmed_comp = trimmed_forest.code
        result[trimmed_comp] = result.get(trimmed_comp, 0) + coeff
        result[comp] = result.get(comp, 0) + coeff * beta * T((-desc,))
    return result

def grove_as_forest_dict(comp, beta, length):
    thmp = grove_as_thompson(comp, beta, length)
    print(f"Thompson result for {comp}: {thmp}")    
    result = {}
    for comp, coeff in thmp.items():
        comp0 = tuple([c for c in comp if c > 0])
        forest_comp = pad_tuple(indexed_forest_from_trimming_word(comp0).code, length)
        result[forest_comp] = result.get(forest_comp, 0) + coeff
    return result


def groth_as_rc(perm, beta, length):
    dct = WCGraph.groth_to_schub(perm, beta=beta)
    result = f.zero
    for perm, coeff in dct.items():
        # for rc in RCGraph.all_rc_graphs(perm, length):
        #     if rc.is_forest_rc:
        result += f.from_dict(dict.fromkeys(RCGraph.all_rc_graphs(perm, length), coeff), snap=True)
    return result

def _forest_weight(wc):
    if wc.perm.inv == len(wc.perm_word):
        return RCGraph(wc).forest_weight
    return wc._snap_reduced().forest_weight

def grove_rc_try(comp, beta, length):
    """GrovePoly polynomial of ``comp`` built by inverting the omega insertion.

    We enumerate the compatible set-valued labelings ``kappa`` of the indexed
    forest ``F = weak_composition_to_indfor(comp)`` (the grove definition), and
    realize each labeling as a WCGraph by writing down its (word, compatible
    sequence) pair *explicitly* -- the published inverse of the set-valued omega
    insertion -- and feeding it to ``WCGraph.from_word_compatible``.

    Concretely:

    * The canonical left-binary-search labeling ``P`` of ``F`` and the decreasing
      labeling ``Q`` of the principal reduced RC graph form the omega pair for
      ``F``. The map ``Gamma = omega_reduced_word_from_labelings`` reads the
      reduced word ``W`` off ``(P, Q)`` -- this is the inverse of the insertion,
      so ``W`` is determined by ``F`` (not chosen at random). Node ``v`` sits at
      word position ``len(W) - Q(v)`` and carries the reduced letter
      ``ell(v) = W[len(W) - Q(v)]``.
    * A labeling ``kappa`` places the letter ``ell(v)`` into every row
      ``r in kappa(v)``. Reading the resulting ``(row, letter)`` pairs in
      ``(row, -letter)`` order gives a weakly-increasing compatible sequence and
      a word whose rows are strictly decreasing, i.e. exactly the data
      ``WCGraph.from_word_compatible`` consumes. The multiset of compatible
      values is the disjoint union of the ``kappa(v)``, so the WCGraph monomial
      equals ``x^kappa``.

    Each WCGraph contributes ``beta**(|kappa| - |F|)`` times its monomial; beta is
    the degree ``-1`` homogenizer recording extra labels beyond one per node.
    """
    from itertools import combinations

    from schubmult.combinatorics.rc_graph import RCGraph

    forest = weak_composition_to_indfor(comp)

    def labelings_below(node, min_value):
        results = []
        for subset_size in range(1, node.rho - min_value + 2):
            for subset in combinations(range(min_value, node.rho + 1), subset_size):
                top = subset[-1]
                left_opts = labelings_below(node.left, top) if node.left is not None else [{}]
                right_opts = labelings_below(node.right, top + 1) if node.right is not None else [{}]
                for left in left_opts:
                    for right in right_opts:
                        labeling = {node.index: tuple(subset)}
                        labeling.update(left)
                        labeling.update(right)
                        results.append(labeling)
        return results

    all_labelings = [{}]
    for root in forest._roots:
        root_labelings = labelings_below(root, 1)
        all_labelings = [{**base, **choice} for base in all_labelings for choice in root_labelings]

    # Omega pair (P, Q) for F, read off the principal reduced RC graph.
    principal = RCGraph.principal_rc(uncode(comp), length)
    p_labeling, q_labeling = principal.omega_invariant

    # Gamma (the published inverse of the insertion) reconstructs the reduced word
    # from (P, Q); reversing matches the forward convention of perm_word.
    word = list(reversed(omega_reduced_word_from_labelings(p_labeling, q_labeling)))

    # The reduced letter carried by node v is W at v's Q-position.
    letter_of = {
        node.index: word[len(word) - q_labeling(node.index)]
        for node in p_labeling.forest.inorder_traversal
    }

    result = 0
    wc_set = set()
    for labeling in all_labelings:
        # Explicit (word, compatible sequence): node v's letter into each row of kappa(v).
        entries = [(row, letter_of[index]) for index, label_set in labeling.items() for row in label_set]
        entries.sort(key=lambda pair: (pair[0], -pair[1]))
        compat_seq = [row for row, _ in entries]
        wc_word = [letter for _, letter in entries]
        wc = WCGraph.from_word_compatible(wc_word, compat_seq, length=length)
        assert wc.is_valid
        wc_set.add(wc)
    for wc in wc_set:
        result += (beta ** (len(wc.perm_word) - sum(comp))) * wc.polyvalue(Sx.genset)
    return result

def _signature(code):
    key = next(iter(br.from_rc_graph(RCGraph.principal_rc(uncode(code)), n).keys()))
    return (len(k) for k in key)

def signature(key):
    return (len(k) for k in key)

def _grove_it_up(comp, bw, n):
    # local imports: keep module-level instances out of the closure that joblib pickles
    from schubmult.rings.polynomial_algebra import GrothendieckPolyBasis, GrovePoly

    grove = 0
    
    # groth = bw.full_groth_elem(uncode(comp), n, 1)
    #grove = 0
    seen = set()
    grippy = GrovePoly(*comp).change_basis(GrothendieckPolyBasis)
    for (groth_perm, _), coeff0 in grippy.items():
        groth = bw.full_groth_elem(groth_perm, n, GrovePoly._basis.beta)
        for key, coeff in groth.items():
            # wc = bw.key_to_wc_graph(key).resize(len(comp))
            # #grove += coeff * bw.full_groth_elem(perm, n, 1)
            # if wc.grove_weight == tuple(comp):# and wc not in seen:
                #assert coeff >= 0
            grove += coeff * coeff0 * bw(key)
            #seen.add(wc)
        # elif wc.perm != uncode(comp):
        #     grove += coeff * bw(key)
        # rc_key = bw.key_to_wc_graph(key).resize(len(comp))
        # if rc_key.grove_weight == tuple(comp) or rc_key.perm != uncode(comp):
        # grove += coeff * bw(key)
    #grove = bw.from_dict({k: v for k, v in grove.items() if (bw.key_to_wc_graph(k).resize(len(comp)).grove_weight == tuple(comp) or bw.key_to_wc_graph(k).perm != uncode(comp))})
    return grove

def _grove_it_up_old(comp, bw, n):
    grove_asgroth = GrovePoly(*comp).change_basis(GrothendieckPolyBasis)
    # groth = bw.full_groth_elem(uncode(comp), n, 1)
    grove = 0
    for (perm, length), coeff in grove_asgroth.items():
        grove += coeff * bw.full_groth_elem(perm, n, 1)
        # rc_key = bw.key_to_wc_graph(key).resize(len(comp))
        # if rc_key.grove_weight == tuple(comp) or rc_key.perm != uncode(comp):
        # grove += coeff * bw(key)
    #grove = bw.from_dict({k: v for k, v in grove.items() if (bw.key_to_wc_graph(k).resize(len(comp)).grove_weight == tuple(comp) or bw.key_to_wc_graph(k).perm != uncode(comp))})
    return grove


_WORKER_STATE = {}


def _worker_state():
    """Per-process ring and grove cache: algebra elements are not picklable, so each
    worker builds its own ``BoundedWCFactorAlgebra`` and memoizes the grove elements it needs."""
    if "bw" not in _WORKER_STATE:
        from schubmult.rings.combinatorial.bounded_wc_factor_algebra import BoundedWCFactorAlgebra

        _WORKER_STATE["bw"] = BoundedWCFactorAlgebra()
        _WORKER_STATE["poles"] = {}
    return _WORKER_STATE["bw"], _WORKER_STATE["poles"]


def _grove_cached(comp, bw, n, poles):
    if comp not in poles:
        poles[comp] = _grove_it_up(comp, bw, n)
    return poles[comp]


def _freeze_elem(elem):
    """Picklable form of a BoundedWCFactorAlgebra element: ``{(factors, size): coeff}``."""
    return {(tuple(key.factors), key.size): coeff for key, coeff in elem.items()}


def _thaw_elem(bw, frozen):
    return bw.from_dict({bw.make_key(factors, size): coeff for (factors, size), coeff in frozen.items()})


def build_grove(comp, n):
    """Phase 1 task: the grove element of ``comp`` in frozen (picklable) form, with its key count."""
    bw, poles = _worker_state()
    elem = _grove_cached(comp, bw, n, poles)
    return comp, _freeze_elem(elem)


def _thawed(comp, frozen, bw, poles):
    if comp not in poles:
        poles[comp] = _thaw_elem(bw, frozen)
    return poles[comp]


def run_single_test(comp1, comp2, n, frozen1=None, frozen2=None):
    """Check the grove LR rule for one pair; returns ``((comp1, comp2), ok, detail)`` (picklable).

    With ``frozen1``/``frozen2`` (from ``build_grove``) the elements are rehydrated instead of
    recomputed; either way they are cached per worker process.
    """
    from schubmult.rings.polynomial_algebra import GrovePoly

    bw, poles = _worker_state()
    grove1 = _thawed(comp1, frozen1, bw, poles) if frozen1 is not None else _grove_cached(comp1, bw, n, poles)
    grove2 = _thawed(comp2, frozen2, bw, poles) if frozen2 is not None else _grove_cached(comp2, bw, n, poles)

    producto = (grove1 * grove2).to_wc_graph_ring_element().resize(n - 1)
    real_prod = GrovePoly(*comp1) * GrovePoly(*comp2)

    checko_prod = 0
    for wc, v in producto.items():
        if wc.grove_weight == wc.length_vector:
            checko_prod += v * GrovePoly(*wc.grove_weight)

    good = real_prod.almosteq(checko_prod)
    detail = None if good else f"{real_prod - checko_prod=}\n{real_prod=}\n{checko_prod=}"
    del producto
    _release_caches(bw)
    return (comp1, comp2), good, detail


def _release_caches(bw):
    """Drop the memo caches that grow with every product, so a worker's footprint is that of its
    largest task rather than the sum of all its tasks.  Everything cleared is a pure memo
    (``WCGraph`` equality is by content, so dropping the interning table is safe too)."""
    from schubmult.combinatorics.wc_graph import WCGraph

    WCGraph.squash_product.cache_clear()
    WCGraph.zero_out_last_row.cache_clear()
    WCGraph.to_mbpd.cache_clear()
    WCGraph.from_mbpd.cache_clear()
    WCGraph.__xnew_cached__.cache_clear()
    WCGraph._z_cache.clear()
    bw._mul_keys.cache_clear()
    bw._normalize_key.cache_clear()


# Memory model for one pair, calibrated at n = 5: the product has about ``|grove1| * |grove2|``
# keys and the squash/MBPD caches hold ~20 KB per key while the task runs (a 47k-key pair peaks
# at 1.1 GB).  ``_release_caches`` drops them after every task, so a worker's footprint is that
# of its heaviest task; chunks group pairs of similar weight and run as many workers as the
# budget allows for the heaviest pair in the chunk.
_GB_PER_KEY = 20e-6
_GB_PER_PROCESS = 0.3


def _pair_gb(sizes, pair):
    return _GB_PER_PROCESS + _GB_PER_KEY * sizes[pair[0]] * sizes[pair[1]]


def _memory_chunks(pairs, sizes, budget_gb, n_jobs):
    """Group ``pairs`` (heaviest first) into chunks; yields ``(chunk, n_jobs_for_chunk)``.

    A chunk's worker count is ``min(n_jobs, budget // heaviest pair)`` for its first (heaviest)
    pair, and it takes at least ``2 * workers`` pairs so no worker idles, continuing while the
    following pairs could not use at least twice that many workers (capped at ``n_jobs``).

    Worker counts are at least 2: ``Parallel(n_jobs=1)`` runs *in the parent*, which would
    both fill the parent with caches and make the task function unpicklable afterwards
    (cloudpickle serializes the module-level worker state by value)."""

    def allowed(pair):
        return max(2, min(n_jobs, int(budget_gb // _pair_gb(sizes, pair))))

    index = 0
    while index < len(pairs):
        jobs = allowed(pairs[index])
        if jobs == n_jobs:
            yield pairs[index:], jobs
            return
        end = index + 1
        while end < len(pairs) and (end - index < 2 * jobs or allowed(pairs[end]) < min(n_jobs, 2 * jobs)):
            end += 1
        yield pairs[index:end], jobs
        index = end


def _available_gb():
    try:
        import psutil

        return psutil.virtual_memory().available / 1024**3
    except ImportError:
        return float(os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_AVPHYS_PAGES")) / 1024**3


if __name__ == "__main__":
    import json
    import os
    import sys
    import time

    from joblib import Parallel, delayed
    from joblib.externals.loky import get_reusable_executor

    n = int(sys.argv[1])
    n_jobs = int(sys.argv[2]) if len(sys.argv) > 2 else os.cpu_count()
    results_path = sys.argv[3] if len(sys.argv) > 3 else f"grove_lr_rule_results_n{n}.jsonl"
    perms = Permutation.all_permutations(n)
    comps = [tuple(perm.pad_code(n -  1)) for perm in perms]

    # Phase 1: every grove element once, in parallel, shipped back in picklable form.
    t0 = time.perf_counter()
    frozen = dict(Parallel(n_jobs=n_jobs, batch_size=1)(delayed(build_grove)(comp, n) for comp in comps))
    sizes = {comp: len(elem) for comp, elem in frozen.items()}
    print(f"Phase 1: {len(comps)} grove elements, {sum(sizes.values())} keys, in {time.perf_counter() - t0:.1f}s "
          f"(largest: {sorted(sizes.items(), key=lambda kv: -kv[1])[:3]})", flush=True)
    get_reusable_executor().shutdown(wait=True)

    # Results are appended per pair so an interrupted run resumes where it stopped.
    done_pairs = {}
    if os.path.exists(results_path):
        with open(results_path) as fh:
            for line in fh:
                rec = json.loads(line)
                done_pairs[(tuple(rec["comp1"]), tuple(rec["comp2"]))] = rec["ok"]
        print(f"Resuming: {len(done_pairs)} pairs already in {results_path}", flush=True)

    # Phase 2: all ordered pairs, largest products first (LPT scheduling), in memory-bounded
    # chunks on fresh workers.  The budget leaves headroom for the interpreters themselves.
    pairs = [pair for pair in sorted(itertools.product(comps, repeat=2), key=lambda pair: -_pair_gb(sizes, pair)) if pair not in done_pairs]
    budget_gb = 0.75 * _available_gb()
    heaviest = _pair_gb(sizes, pairs[0]) if pairs else 0.0
    print(f"Phase 2: {len(pairs)} pairs, memory budget {budget_gb:.1f} GB, heaviest pair ~{heaviest:.1f} GB", flush=True)
    if 2 * heaviest > budget_gb:
        print("WARNING: two of the heaviest pairs exceed the memory budget together; the first chunk may swap", flush=True)

    t0 = time.perf_counter()
    failures = []
    done = 0
    total = len(pairs)
    with open(results_path, "a") as out:
        for chunk, chunk_jobs in _memory_chunks(pairs, sizes, budget_gb, n_jobs):
            print(f"-- chunk of {len(chunk)} pairs on {chunk_jobs} workers (heaviest ~{_pair_gb(sizes, chunk[0]):.1f} GB each)", flush=True)
            for (comp1, comp2), good, detail in Parallel(n_jobs=chunk_jobs, batch_size=1, return_as="generator_unordered")(
                delayed(run_single_test)(comp1, comp2, n, frozen[comp1], frozen[comp2]) for comp1, comp2 in chunk
            ):
                done += 1
                out.write(json.dumps({"comp1": comp1, "comp2": comp2, "ok": bool(good)}) + "\n")
                out.flush()
                if good:
                    print(f"Success {comp1} {comp2}  [{done}/{total}, {time.perf_counter() - t0:.0f}s]", flush=True)
                else:
                    failures.append(((comp1, comp2), detail))
                    print(f"Failed for {comp1} * {comp2}: {detail}", flush=True)
            # fresh worker processes for the next chunk: drops every per-process cache
            get_reusable_executor().shutdown(wait=True)

    previous_failures = [pair for pair, ok in done_pairs.items() if not ok]
    print(f"{total - len(failures)}/{total} pairs succeeded in {time.perf_counter() - t0:.0f}s "
          f"({len(done_pairs)} from a previous run, {len(previous_failures)} of them failed)", flush=True)
    if failures or previous_failures:
        print("Failures:", [key for key, _ in failures] + previous_failures)
    sys.exit(1 if failures or previous_failures else 0)
