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
    from schubmult import BoundedRCFactorAlgebra, RCGraph, RCGraphRing, Sx, uncode, CrystalGraphTensor
    from schubmult.utils.tuple_utils import pad_tuple
    from sympy import Add, Mul, expand, Pow, sympify, pretty_print

    g = BoundedRCFactorAlgebra()
    r = RCGraphRing()

    # def cem_schub(perm, n):
    #     return sum([g.from_tensor_dict(cem_dict, n) for rc, cem_dict in RCGraph.full_CEM(perm, n).items()])

    # def cem_schub_schur_decomp(perm, n):
        
    #     result = Sx.zero @ Sx.zero
    #     cd = ((perm.strict_mul_dominant(n))).trimcode
    #     if any(a < n for a in cd):
    #         toadd = min(n - a for a in cd if a < n)
    #         cd = [a + toadd for a in cd]
    #     domperm = uncode(cd)
    #     reppy = sympify(expand(Sx(perm).cem_rep(mumu=~domperm, elem_func=Sx.symbol_elem_func), func=False))
    #     for arg in Add.make_args(reppy):
    #         coeff, schur_part = arg.as_coeff_Mul()
    #         part1 = Sx.one
    #         part2 = Sx.one
    #         for elem_arg in Mul.make_args(schur_part):
    #             if isinstance(elem_arg, Pow):
    #                 base, exp = elem_arg.as_base_exp()
    #             else:
    #                 base = elem_arg
    #                 exp = 1
    #             for _ in range(exp):
    #                 if base.numvars < n:
    #                     part1 *= base
    #                 else:
    #                     part2 *= base
    #             # if elem_arg.numvars < n:
    #             #     part1 *= elem_arg
    #             # else:
    #             #     part2 *= elem_arg
    #         result += coeff * part1 @ part2
    #     return result

    try:
        # g_result = g.zero
        result = r.zero
        
        length = max(len(perm1), len(perm2))
        # length = n
        #max(len(perm1.trimcode), len(perm2.trimcode)) + 1
        #prd = Sx.from_dict({k: v for k, v in (Sx(perm1) * Sx(perm2)).items() if len(k) <= length})
        prd = Sx(perm1) * Sx(perm2)
        #length = n
        # partition1 = tuple((~(perm1.strict_mul_dominant(length))).trimcode)
        # partition2 = tuple((~(perm2.strict_mul_dominant(length))).trimcode)
        # decomp1 = cem_schub_schur_decomp(perm1, length)
        schub1 = g.full_schub_elem(perm1, length)
        schub2 = g.full_schub_elem(perm2, length)
        
        rc_by_grass1 = (g@g).zero
        rc_by_grass2 = (g@g).zero
        for key1, coeff1 in schub1.items():
            if len(key1) > 0 and len(key1[-1]) == key1.size:
                rc_by_grass1 += coeff1 * g(g.make_key(key1[:-1], key1.size)) @ g(g.make_key((key1[-1],), key1.size))
            else:
                rc_by_grass1 += coeff1 * g(g.make_key(key1, key1.size)) @ g(g.make_key((), key1.size))

        for key2, coeff2 in schub2.items():
            if len(key2) > 0 and len(key2[-1]) == key2.size:
                rc_by_grass2 += coeff2 * g(g.make_key(key2[:-1], key2.size)) @ g(g.make_key((key2[-1],), key2.size))
                #continue
            else:
                rc_by_grass2 += coeff2 * g(g.make_key(key2, key2.size)) @ g(g.make_key((), key2.size))
        overflow_result = (g@g).zero

        for (part1, grass1), coeff1 in rc_by_grass1.items():
            for (part2, grass2), coeff2 in rc_by_grass2.items():
                base = g(part1) * g(part2)
                grass = g(grass1) * g(grass2)
                overflow_result += base @ grass
                # for key_base, coeff in base.items():
                #     if len(key_base) > 0 and len(key_base[-1]) == key_base.size:
                #         overflow_result += coeff * coeff1 * coeff2 * g(g.make_key(key_base[:-1], key_base.size)) @ (g(g.make_key((key_base[-1],), key_base.size)) * g(grass1) * g(grass2))
                #         #continue
                #     else:
                #         overflow_result += coeff * coeff1 * coeff2 * g(g.make_key(key_base, key_base.size)) @ (g(grass1) * g(grass2))
                #overflow_result += coeff1 * coeff2 * g(g.make_key((grass1,), length)) @ g(g.make_key((grass2,), length)) @ (g(part1)*g(part2)).to_rc_graph_ring_element()
        #pretty_print(overflow_result)
        # for key2, coeff2 in schub2.items():
        # # print(f"Grassmann decomposition for {perm1}:")
        # # for grass, part in rc_by_grass1.items():
        # #     print(f"  Grassmann part {grass}: {part}")
        # # print(f"Grassmann decomposition for {perm2}:")
        # # for grass, part in rc_by_grass2.items():
        # #     print(f"  Grassmann part {grass}: {part}")
        # #input()
        # #return True
        # rc_by_grass_result = {}
        # tensor_result = (r@g@g).zero
        # for grass1, part1 in rc_by_grass1.items():
        #     for grass2, part2 in rc_by_grass2.items():
        #         tensor_result += (g(g.make_key((grass1,), length))* g(g.make_key((grass2,), length))).to_rc_graph_ring_element() @ part1 @ part2
        #         # part_result = (part1 * part2)
        # tensor_result2 = (r@r).zero 
        # for (grass, part1, part2), coeff in tensor_result.items():
        #     tensor_result2 += coeff * r(grass) @ (g(part1)*g(part2)).to_rc_graph_ring_element()
        #         # if part_result != g.zero:
        #             #combined_grass = grass1.squash_product(grass2)
        #             #rc_by_grass_result[combined_grass] = rc_by_grass_result.get(combined_grass, 0) + part_result
        #             # for key2, coeff2 in part_result.items():
        #             #     key_grass = combined_grass
        #             #     if len(key2) > 0 and len(key2[-1]) == key2.size:
        #             #         key_grass = key_grass.squash_product(key2[-1])
                            
        #             #         rc_by_grass_result[key_grass] = rc_by_grass_result.get(key_grass, 0) + coeff2 * g(g.make_key(key2[:-1], key2.size))
        #             #     else:
        #             #         rc_by_grass_result[key_grass] = rc_by_grass_result.get(key_grass, 0) + coeff2 * g(g.make_key(key2, key2.size))
        result = r.zero
        
        for (part, grass), coeff in overflow_result.items():
            # if grass.perm.inv != 0:
            #     continue
            #assert coeff >= 0
            result += coeff * (g(part)*g(grass)).to_rc_graph_ring_element()#.resize(len(grass)).squash_product(grass))
        
        # for grass, part in rc_by_grass_result.items():
        #     # if not grass.is_highest_weight:
        #     #     continue
        #     # hw_sum = g.zero
        #     # for part2, coeff in part.items():
        #     #     if part2.is_highest_weight:
        #     #         if CrystalGraphTensor(*part2, grass).is_highest_weight:
        #     #         # print("Grass")
        #     #         # pretty_print(grass)
        #     #         # print("Part")
        #     #         # pretty_print(part)
        #     #             hw_sum += coeff * g(part2)
        #     # print("Grass")
        #     # pretty_print(grass)
        #     # print("Part")
        #     # pretty_print(hw_sum)
        #     result += part * g(g.make_key((grass,), length))
        # result = result.to_rc_graph_ring_element()
        # # pretty_print(result.to_rc_graph_ring_element().resize(length))
        # # return True
        # # result = g_result.to_rc_graph_ring_element().resize(n)
        # # if any(v < 0 for v in result.values()):
        # #     print(f"Negative coefficient in result for {perm1} and {perm2}: {result}")
        # #     return False
        
        # pretty_print(result)
        prd2 = Sx.zero
        for rc, coeff in result.items():
            if coeff != prd.get(rc.perm, 0):
                print(f"Coeff mismatch for {perm1}, {perm2} at {rc.perm}: got {coeff}, expected {prd.get(rc.perm, 0)}")
                return False
            if rc.is_principal:#is_highest_weight and rc.extremal_weight == pad_tuple(rc.perm.trimcode, len(rc)):
                # if rc.perm not in seen:
                #     seen.add(rc.perm)
                prd2 += coeff * Sx(rc.perm)

        if prd != prd2:
            print(f"Product mismatch for {perm1}, {perm2}: expected {prd}, got {prd2}")
            return False
        #from sympy import pretty_print
        #pretty_print(result)
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
        if skip_id and (perm2.inv != 0):# or perm2.inv == 0):
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

