import argparse
import gc
import itertools
import json
import os
import shutil
import time
import uuid
import inspect
import sys
from joblib import Parallel, delayed


def _as_jsonable_key_part(value):
    """Convert key parts into JSON-safe values.

    Lists/tuples are serialized as lists so tuple-based keys survive JSON round-trips.
    """
    if isinstance(value, (tuple, list)):
        return [_as_jsonable_key_part(v) for v in value]
    return value


def _as_hashable_key_part(value):
    """Convert JSON-loaded key parts back into hashable Python values."""
    if isinstance(value, list):
        return tuple(_as_hashable_key_part(v) for v in value)
    return value


def _key_to_str(k):
    """Tuple key -> JSON-safe string key.

    Supports tuple keys of any arity.
    """
    return json.dumps([_as_jsonable_key_part(part) for part in k])


def _str_to_key(s):
    """Inverse of _key_to_str -> hashable tuple key of any supported arity."""
    parts = json.loads(s)
    return tuple(_as_hashable_key_part(part) for part in parts)


def safe_save_recording(obj, filename, meta=None):
    temp_json = f"{filename}.json.tmp"
    json_file = f"{filename}.json"
    try:
        payload = {"meta": meta or {}, "records": {_key_to_str(k): v for k, v in obj.items()}}
        with open(temp_json, "w") as f:
            json.dump(payload, f)
        if os.path.exists(json_file):
            shutil.copy2(json_file, f"{json_file}.backup")
        os.replace(temp_json, json_file)
    except Exception:
        import traceback
        traceback.print_exc()
        if os.path.exists(temp_json):
            os.remove(temp_json)


def randomized_backup_copy(path, tag="backup"):
    if not os.path.exists(path):
        return None
    stamp = time.strftime("%Y%m%d-%H%M%S")
    backup_path = f"{path}.{tag}.{stamp}.{uuid.uuid4().hex[:8]}"
    shutil.copy2(path, backup_path)
    return backup_path


def safe_load_recording(filename):
    json_file = f"{filename}.json"
    if os.path.exists(json_file):
        try:
            with open(json_file) as f:
                loaded = json.load(f)
            if isinstance(loaded, dict) and "records" in loaded:
                raw_records = loaded.get("records", {})
                meta = loaded.get("meta", {}) or {}
            else:
                raw_records = loaded
                meta = {}
            dct = {}
            for k, v in raw_records.items():
                tp = _str_to_key(k)
                dct[tp] = v
            return dct, meta
        except Exception as exc:
            backup_path = randomized_backup_copy(json_file, tag="parsefail")
            print(f"Failed to parse existing verification JSON: {json_file}", flush=True)
            print(f"Parse error: {exc}", flush=True)
            if backup_path:
                print(f"Copied unreadable file to: {backup_path}", flush=True)
            return {}, {"parse_failed": True, "parse_backup": backup_path}
    return {}, {}


def verify_pair(perm1, perm2, n):
    from schubmult import Sx, uncode
    from schubmult.rings.combinatorial.forest_rc_ring import ForestRCGraphRing
    from schubmult.utils.tuple_utils import pad_tuple
    from sympy import Add, Mul, expand, Pow, sympify, pretty_print
    from schubmult.rings.polynomial_algebra import PolynomialAlgebra, ForestPolyBasis

    ForestPoly = PolynomialAlgebra(ForestPolyBasis(Sx.genset))

    
    r = ForestRCGraphRing()

    

    try:
        result = r.zero
        
        prd = 0
    
        length = max(len(perm1), len(perm2)) - 1
    
        prd = ForestPoly(perm1.pad_code(length)) * ForestPoly(perm2.pad_code(length))
        
        comp1 = perm1.pad_code(length)
        comp2 = perm2.pad_code(length)
        forest1 = r.forest_poly(comp1)
        forest2 = r.forest_poly(comp2)
        if any(v < 0 for v in forest1.values()):
            print(f"Negative coefficient in forest1 for {perm1}: {forest1}")
            return False
        if any(v < 0 for v in forest2.values()):
            print(f"Negative coefficient in forest2 for {perm2}: {forest2}")
            return False
        result = r.dual_product(forest1, forest2)

        if any(v < 0 for v in result.values()):
            print(f"Negative coefficient in result for {perm1} and {perm2}: {result}")
            return False

        prd2 = 0
        for rc, coeff in result.items():
            if coeff != prd.get(rc.forest_weight, 0):
                print(f"Coeff mismatch for {perm1}, {perm2} at {rc.forest_weight}: got {coeff}, expected {prd.get(rc.forest_weight, 0)}")
                return False
            
            if rc.forest_weight == rc.length_vector:
                prd2 += coeff * ForestPoly(rc.forest_weight)

        
        if not prd.almosteq(prd2):
            print(f"Product mismatch for {perm1}, {perm2}: expected {prd}, got {prd2}")
            return False

        return True
    except Exception as e:
        import traceback
        print(f"Exception verifying ({perm1}, {perm2}, {n}): {traceback.format_exc()}")
        return False


def clear_rcgraph_caches():
    """Clear all known RCGraph and related caches to reduce memory overhead."""
    from schubmult.combinatorics import rc_graph
    RCGraph = rc_graph.RCGraph
    cleared = set()
    # Clear all cache-decorated methods
    for name, obj in inspect.getmembers(RCGraph):
        if hasattr(obj, "cache_clear") and callable(obj.cache_clear):
            try:
                obj.cache_clear()
                cleared.add(name)
            except Exception:
                pass
    # Explicitly clear static-level caches
    static_caches = [
        "_graph_cache",
        "_cache_by_weight",
        "_z_cache",
        "w_key_cache",
        "rc_cache",
        "_rc_cache",
    ]
    for cache_name in static_caches:
        if hasattr(RCGraph, cache_name):
            try:
                cache_obj = getattr(RCGraph, cache_name)
                if hasattr(cache_obj, "clear"):
                    cache_obj.clear()
                elif isinstance(cache_obj, set):
                    cache_obj.clear()
                elif isinstance(cache_obj, dict):
                    cache_obj.clear()
                cleared.add(cache_name)
            except Exception:
                pass
    # Clear module-level caches if any
    for name, obj in inspect.getmembers(sys.modules[rc_graph.__name__]):
        if hasattr(obj, "cache_clear") and callable(obj.cache_clear):
            try:
                obj.cache_clear()
                cleared.add(name)
            except Exception:
                pass
    return cleared


def run_single_test(perm1, perm2, n):
    key = (tuple(perm1), tuple(perm2), n)
    good = verify_pair(perm1, perm2, n)
    status = "Success" if good else "FAIL"
    print(f"{status} ({perm1.trimcode}, {perm2.trimcode}, {n}) at {time.ctime()}", flush=True)
    return key, good


def recording_saver(shared_updates_dict, lock, verification_filename, stop_event, base_recording=None, meta=None, sleep_time=10):
    full_recording = dict(base_recording or {})
    while not stop_event.is_set():
        with lock:
            keys = list(shared_updates_dict.keys())
            for k in keys:
                try:
                    full_recording[k] = shared_updates_dict[k]
                except Exception:
                    continue
            for k in keys:
                try:
                    del shared_updates_dict[k]
                except Exception:
                    continue

        if keys:
            print(f"Saving {len(full_recording)} entries to {verification_filename} at {time.ctime()} ({len(keys)} new)", flush=True)
            safe_save_recording(full_recording, verification_filename, meta=meta or {})
            gc.collect()

        time.sleep(sleep_time)

    # Final flush and save
    with lock:
        keys = list(shared_updates_dict.keys())
        for k in keys:
            try:
                full_recording[k] = shared_updates_dict[k]
            except Exception:
                continue
        for k in keys:
            try:
                del shared_updates_dict[k]
            except Exception:
                continue

    safe_save_recording(full_recording, verification_filename, meta=meta or {})
    gc.collect()



def generate_test_cases(perms, n, skip_id, completed_true_keys):
    left_factor = perms
    right_factor = perms
    product_set = set(itertools.product(left_factor, right_factor))
    product_set = {key for key in product_set if (tuple(key[0]), tuple(key[1]), n) not in completed_true_keys}
    for perm1, perm2 in product_set:
        if skip_id and (perm1.inv == 0 or perm2.inv == 0):
            continue
        key = (tuple(perm1), tuple(perm2), n)
        if key in completed_true_keys:
            continue
        yield perm1, perm2, n


def _merge_with_updates(base_recording, shared_updates_dict, lock):
    merged = dict(base_recording)
    with lock:
        keys = list(shared_updates_dict.keys())
    for k in keys:
        try:
            merged[k] = shared_updates_dict[k]
        except Exception:
            continue
    return merged


def main():
    from schubmult import Permutation

    parser = argparse.ArgumentParser(description="Final LR rule verification (multiprocessed)")
    parser.add_argument("n", type=int, help="Permutation size parameter")
    parser.add_argument("num_processors", type=int, help="Number of worker processes")
    parser.add_argument("filename", nargs="?", help="Base filename for verification output (json will be <filename>.verification.json)")
    parser.add_argument("--skip-id", action="store_true", default=True, help="Skip identity permutations (default: True)")
    parser.add_argument("--no-skip-id", dest="skip_id", action="store_false", help="Include identity permutations")

    args = parser.parse_args()
    n = args.n
    num_processors = args.num_processors
    skip_id = args.skip_id
    filename = args.filename
    verification_filename = filename + ".verification" if filename else None

    meta = {"n": n, "skip_id": skip_id}

    loaded_recording = {}
    load_meta = {}
    if verification_filename:
        json_file = f"{verification_filename}.json"
        if os.path.exists(json_file):
            preload_backup = randomized_backup_copy(json_file, tag="preload")
            if preload_backup:
                print(f"Pre-load backup created: {preload_backup}", flush=True)
        loaded_recording, load_meta = safe_load_recording(verification_filename)
        if load_meta.get("parse_failed"):
            print("Parse failed for existing verification file; continuing without persistence to avoid overwrite.", flush=True)
            verification_filename = None
            loaded_recording = {}

    perms = Permutation.all_permutations(n)
    #perms.sort(key=lambda p: (p.inv, p.trimcode))

    print(f"Final LR rule verification: n={n}, {len(perms)} permutations, {num_processors} workers", flush=True)

    completed_true_keys = {k for k, v in loaded_recording.items() if v is True}


    if loaded_recording:
        already_done = sum(1 for v in loaded_recording.values() if v is True)
        print(f"Loaded {len(loaded_recording)} existing results ({already_done} verified)", flush=True)

    # Prepare test cases
    test_cases = list(generate_test_cases(perms, n, skip_id, completed_true_keys))
    total = len(test_cases)
    print(f"Running {total} test cases using joblib with {num_processors} workers", flush=True)

    # Run in parallel with joblib, saving after each batch
    batch_size = 10
    final_recording = dict(loaded_recording)
    total = len(test_cases)
    for batch_start in range(0, total, batch_size):
        batch = test_cases[batch_start:batch_start+batch_size]
        results = Parallel(n_jobs=num_processors, batch_size=1, verbose=10)(
            delayed(run_single_test)(perm1, perm2, n) for (perm1, perm2, n) in batch
        )
        for key, good in results:
            final_recording[key] = good
        if verification_filename:
            safe_save_recording(final_recording, verification_filename, meta=meta)
        print(f"Saved progress: {min(batch_start+batch_size, total)}/{total} cases complete", flush=True)

    if verification_filename:
        final_recording, _ = safe_load_recording(verification_filename)

    # Report
    total = len(final_recording)
    failures = sum(1 for v in final_recording.values() if v is False)
    successes = sum(1 for v in final_recording.values() if v is True)

    if failures:
        print(f"\n{failures} FAILURES out of {total}:")
        for k, v in final_recording.items():
            if v is False:
                print(f"  {k}")
    else:
        print(f"\nAll {successes} pairs verified successfully!")


if __name__ == "__main__":
    main()

