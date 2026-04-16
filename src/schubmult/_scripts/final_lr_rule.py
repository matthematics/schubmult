import argparse
import gc
import itertools
import json
import os
import shutil
import time
from multiprocessing import Event, Manager, Process


def _key_to_str(k):
    """(perm1, perm2, n) -> JSON-safe string key using list representations."""
    perm1, perm2, n = k
    return json.dumps([list(perm1), list(perm2), n])


def _str_to_key(s):
    """Inverse of _key_to_str -> tuple of (tuple, tuple, int)."""
    perm1, perm2, n = json.loads(s)
    return (tuple(perm1), tuple(perm2), n)


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


def safe_load_recording(filename):
    json_file = f"{filename}.json"
    if os.path.exists(json_file):
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
    return {}, {}


def verify_pair(perm1, perm2, n):
    from schubmult import BoundedRCFactorAlgebra, RCGraph, RCGraphRing, Sx, uncode
    from schubmult.utils.tuple_utils import pad_tuple
    from sympy import Add, Mul, expand

    g = BoundedRCFactorAlgebra()
    r = RCGraphRing()

    # def cem_schub(perm, n):
    #     return sum([g.from_tensor_dict(cem_dict, n) for rc, cem_dict in RCGraph.full_CEM(perm, n).items()])

    # def cem_schub_schur_decomp(perm, n):
    #     result = Sx.zero @ Sx.zero
    #     cd = perm.strict_mul_dominant().trimcode
    #     if any(a < n for a in cd):
    #         toadd = min(n - a for a in cd if a < n)
    #         cd = [a + toadd for a in cd]
    #     domperm = uncode(cd)
    #     reppy = expand(Sx(perm).cem_rep(mumu=~domperm, elem_func=Sx.symbol_elem_func), func=False)
    #     for arg in Add.make_args(reppy):
    #         coeff, schur_part = arg.as_coeff_Mul()
    #         part1 = Sx.one
    #         part2 = Sx.one
    #         for elem_arg in Mul.make_args(schur_part):
    #             if elem_arg.numvars < n:
    #                 part1 *= elem_arg
    #             else:
    #                 part2 *= elem_arg
    #         result += coeff * part1 @ part2
    #     return result

    try:
        g_result = g.zero
        result = r.zero
        prd = Sx(perm1) * Sx(perm2)
        length = max(len(perm1), len(perm2)) - 1
        partition1 = tuple((~(perm1.strict_mul_dominant(length))).trimcode)
        partition2 = tuple((~(perm2.strict_mul_dominant(length))).trimcode)
        schub1 = g.schub_elem(perm1, length, partition=partition1)
        #schub1 = g.from_tensor_dict(schub1_base, size=length)
        schub2 = g.schub_elem(perm2, length, partition=partition2)
        #schub2 = g.from_tensor_dict(schub2_base, size=length)
        tensor_result = r.zero @ r.zero
        for key1, coeff1 in schub1.items():
            # rc1 = next(iter(g(key1).to_rc_graph_ring_element().resize(n)))
            # if rc1.perm != perm1:
            #     continue
            sumup = r.zero
            # if not key1.is_highest_weight:
            #     continue
            for key2, coeff2 in schub2.items():
                # rc2 = next(iter(g(key2).to_rc_graph_ring_element().resize(n)))
                # if rc2.perm != perm2:
                #     continue
                graph_base = (g(key1) * g(key2))
                for the_key, _ in graph_base.items():
                    if the_key.is_highest_weight:
                        rc = next(iter(g(the_key).to_rc_graph_ring_element().resize(n)))
                        # if prd.get(rc.perm, 0) != 0:
                        #     # if coeff1 * coeff2 < 0:
                        #     #     print(f"Negative coefficient for {rc} in product of {perm1} and {perm2} with keys {key1} and {key2}")
                        #     #     return False
                        #     tensor_result += coeff1 * coeff2 * g(key1).to_rc_graph_ring_element() @ g(key2).to_rc_graph_ring_element()
                        # g_result += coeff1 * coeff2 * g(the_key)
                        sumup += coeff1 * coeff2 * r(rc)
            # if any(v < 0 for v in sumup.values()):
            #     print(f"Negative coefficient in intermediate sumup for {perm1} and {perm2} at key {key1}: {sumup}")
            #     return False
            result += sumup
                        

        # if any(v < 0 for v in result.values()):
        #     print(f"Negative coefficient in result for {perm1} and {perm2}: {result}")
        #     return False
        # if any(v < 0 for v in g_result.values()):
        #     print(f"Negative coefficient in result for {perm1} and {perm2}: {g_result}")
        #     return False
        # result = g_result.to_rc_graph_ring_element().resize(n)
        prd2 = Sx.zero
        for rc, coeff in result.items():
            if coeff != prd.get(rc.perm, 0):
                print(f"Coeff mismatch for {perm1}, {perm2} at {rc.perm}: got {coeff}, expected {prd.get(rc.perm, 0)}")
                return False
            if rc.extremal_weight == pad_tuple(rc.perm.trimcode, len(rc)):
                prd2 += coeff * Sx(rc.perm)

        if prd != prd2:
            print(f"Product mismatch for {perm1}, {perm2}: expected {prd}, got {prd2}")
            return False

        return True
    except Exception as e:
        import traceback
        print(f"Exception verifying ({perm1}, {perm2}, {n}): {traceback.format_exc()}")
        return False


def worker(shared_recording_dict, lock, task_queue):
    while True:
        try:
            item = task_queue.get()
            if item is None:
                break
            perm1, perm2, n = item
        except Exception:
            import traceback
            traceback.print_exc()
            break

        key = (tuple(perm1), tuple(perm2), n)
        with lock:
            if shared_recording_dict.get(key) is True:
                continue

        good = verify_pair(perm1, perm2, n)

        with lock:
            shared_recording_dict[key] = good

        status = "Success" if good else "FAIL"
        print(f"{status} ({perm1.trimcode}, {perm2.trimcode}, {n}) at {time.ctime()}", flush=True)


def recording_saver(shared_recording_dict, lock, verification_filename, stop_event, meta=None, sleep_time=10):
    last_len = -1
    while not stop_event.is_set():
        try:
            new_len = len(shared_recording_dict)
        except Exception:
            new_len = last_len

        if new_len > last_len:
            last_len = new_len
            print(f"Saving {new_len} entries to {verification_filename} at {time.ctime()}", flush=True)
            with lock:
                keys = list(shared_recording_dict.keys())
            recording_copy = {}
            for k in keys:
                try:
                    recording_copy[k] = shared_recording_dict[k]
                except Exception:
                    continue
            safe_save_recording(recording_copy, verification_filename, meta=meta or {})

        time.sleep(sleep_time)

    # Final save
    with lock:
        keys = list(shared_recording_dict.keys())
    recording_copy = {}
    for k in keys:
        try:
            recording_copy[k] = shared_recording_dict[k]
        except Exception:
            continue
    safe_save_recording(recording_copy, verification_filename, meta=meta or {})


def queue_producer(task_queue, perms, n, num_processors, skip_id, grass_left_factor, grass_right_factor, shared_recording_dict):
    task_count = 0
    left_factor = perms
    right_factor = perms
    bad_patterns = [[4,1,3,2],[1,4,3,2],[3,1,4,2]]
    if grass_left_factor:
        left_factor = [p for p in perms if all(not p.has_pattern(pat) for pat in bad_patterns)]
    if grass_right_factor:
        right_factor = [p for p in perms if all(not p.has_pattern(pat) for pat in bad_patterns)]
    for perm1, perm2 in itertools.product(left_factor, right_factor):
        if skip_id and (perm1.inv == 0 or perm2.inv == 0):
            continue
        key = (tuple(perm1), tuple(perm2), n)
        if shared_recording_dict.get(key) is True:
            continue
        task_queue.put((perm1, perm2, n))
        task_count += 1
        if task_count % 100 == 0:
            gc.collect()

    for _ in range(num_processors):
        task_queue.put(None)


def main():
    from schubmult import Permutation

    parser = argparse.ArgumentParser(description="Final LR rule verification (multiprocessed)")
    parser.add_argument("n", type=int, help="Permutation size parameter")
    parser.add_argument("num_processors", type=int, help="Number of worker processes")
    parser.add_argument("filename", nargs="?", help="Base filename for verification output (json will be <filename>.verification.json)")
    parser.add_argument("--skip-id", action="store_true", default=True, help="Skip identity permutations (default: True)")
    parser.add_argument("--grass_left_factor", action="store_true", default=False, help="Use Grassmann left factor (default: False)")
    parser.add_argument("--grass_right_factor", action="store_true", default=False, help="Use Grassmann right factor (default: False)")
    parser.add_argument("--no-skip-id", dest="skip_id", action="store_false", help="Include identity permutations")
    parser.add_argument("--extra", type=int, default=1, help="Extra size for permutation generation (default: 1)")

    args = parser.parse_args()
    n = args.n
    num_processors = args.num_processors
    extra = args.extra
    skip_id = args.skip_id
    grass_left_factor = args.grass_left_factor
    grass_right_factor = args.grass_right_factor
    filename = args.filename
    verification_filename = filename + ".verification" if filename else None

    meta = {"n": n, "extra": extra, "skip_id": skip_id}

    loaded_recording = {}
    if verification_filename:
        loaded_recording, _ = safe_load_recording(verification_filename)

    perms = [p for p in Permutation.all_permutations(n + extra) if len(p.trimcode) <= n]
    perms.sort(key=lambda p: (p.inv, p.trimcode))

    print(f"Final LR rule verification: n={n}, extra={extra}, {len(perms)} permutations, {num_processors} workers", flush=True)

    with Manager() as manager:
        shared_recording_dict = manager.dict()
        lock = manager.Lock()
        stop_event = Event()

        if loaded_recording:
            shared_recording_dict.update(loaded_recording)
            already_done = sum(1 for v in loaded_recording.values() if v is True)
            print(f"Loaded {len(loaded_recording)} existing results ({already_done} verified)", flush=True)

        if verification_filename:
            saver_proc = Process(
                target=recording_saver,
                args=(shared_recording_dict, lock, verification_filename, stop_event, meta),
            )
            saver_proc.start()

        task_queue = manager.Queue()

        producer_proc = Process(
            target=queue_producer,
            args=(task_queue, perms, n, num_processors, skip_id, grass_left_factor, grass_right_factor, shared_recording_dict),
        )
        producer_proc.start()

        workers = []
        for _ in range(num_processors):
            p = Process(target=worker, args=(shared_recording_dict, lock, task_queue))
            p.start()
            workers.append(p)

        producer_proc.join()
        for p in workers:
            p.join()

        stop_event.set()
        if verification_filename:
            saver_proc.join()

        # Report
        total = len(shared_recording_dict)
        failures = sum(1 for v in shared_recording_dict.values() if v is False)
        successes = sum(1 for v in shared_recording_dict.values() if v is True)

        if failures:
            print(f"\n{failures} FAILURES out of {total}:")
            for k, v in shared_recording_dict.items():
                if v is False:
                    print(f"  {k}")
        else:
            print(f"\nAll {successes} pairs verified successfully!")


if __name__ == "__main__":
    main()

