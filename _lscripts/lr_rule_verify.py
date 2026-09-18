import argparse
import base64
import pickle
import shutil
import sys
import time
import os
import json
from multiprocessing import Event, Manager, Process
import gc


def json_key_rep(k):
    """
    Serialize an arbitrary Python object `k` to a JSON-key-safe string.
    Uses pickle + base64 with a 'PICKLE:' prefix; falls back to repr.
    """
    try:
        raw = pickle.dumps(k, protocol=pickle.HIGHEST_PROTOCOL)
        return "PICKLE:" + base64.b64encode(raw).decode("ascii")
    except Exception:
        return "REPR:" + repr(k)


def json_key_load(s):
    """
    Inverse of json_key_rep.
    """
    if not isinstance(s, str):
        raise TypeError("json_key_load expects a string")
    if s.startswith("PICKLE:"):
        payload = s[len("PICKLE:") :]
        raw = base64.b64decode(payload.encode("ascii"))
        return pickle.loads(raw)
    if s.startswith("REPR:"):
        rep = s[len("REPR:") :]
        return eval(rep)
    try:
        return eval(s)
    except Exception:
        return s


def safe_save_recording(obj, filename, meta=None):
    """
    Save {meta: {...}, records: {json_key: value}} atomically.
    """
    temp_json = f"{filename}.json.tmp"
    json_file = f"{filename}.json"
    try:
        payload = {"meta": meta or {}, "records": {json_key_rep(k): v for k, v in obj.items()}}
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
    """
    Load the verification JSON.
    Returns (records_dict, meta_dict).
    Backwards compatible with old flat mapping format (treats file as records and meta={}).
    """
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
            tp = json_key_load(k)
            dct[tp] = v
        return dct, meta
    return {}, {}


def worker(nn, shared_recording_dict, lock, task_queue, y, z):  # noqa: ARG001
    import sympy

    from schubmult import CrystalGraphTensor, DSx, FreeAlgebra, Permutation, RCGraph, RCGraphRing, RootTableau, SchubertBasis, Sx  # noqa: F401

    
    while True:
        from sympy import S

        try:
            item = task_queue.get()
            if item is None:
                print("Got None")
                break
            (key, check_val) = item
            (u, v, w, n) = key
            ring = Sx
            if y is not None and z is not None:
                ring = DSx([]).ring

            check_val = ring.from_dict(check_val)
        except Exception:
            import traceback
            traceback.print_exc()
            break  # queue empty, exit
        with lock:
            if key in shared_recording_dict:
                if shared_recording_dict[key] is True:
                    # print(f"{(key, check_val)} already verified, returning.")
                    continue
                print(f"Previous failure on {key}, will retry.")

        good = False
        # if True:
        # hw_tab = RCGraph.principal_rc(u, n - 1).to_highest_weight()[0]

    
        sm0 = S.Zero
        dualps = {}
        # for w in check_elem:
            # if w != u:
            #     continue
        #for v_rc in RCGraph.all_rc_graphs(v, n):

        hw_tensors = set()
        dom_rc = RCGraph.principal_rc(u, n)
        for v_rc in RCGraph.all_rc_graphs(v, n):
            tensor = CrystalGraphTensor(dom_rc, v_rc)
            hw_tensors.add(tensor.to_highest_weight()[0])

        visited = set()
        #check_elem = Sx(u) * Sx(v)
        for w in check_val:
            # if not good:
            #     break
            for tensor_hw in hw_tensors:
                # if not good:
                #     break
                for inner_tensor in tensor_hw.full_crystal:
                    # if not good:
                    #     break
                    v_rc = inner_tensor.factors[1]
                    dualpocket = v_rc.dualpieri(u, w)
                    if len(dualpocket) > 0:
                        # if tensor_hw in visited:
                        #     good = False
                        #     break
                        # if len() == 0:
                        #     good = False
                        #     break
                        # rc_high = RCGraph.all_rc_graphs(w, n, weight=tensor_hw.crystal_weight)
                        # rc_low = RCGraph.all_rc_graphs(w, n, weight=inner_tensor.to_lowest_weight()[0].crystal_weight)
                        # good = False
                        # for rc_lw in rc_low:
                        #     if not rc_lw.is_lowest_weight or rc_lw.to_highest_weight()[0] not in rc_high:
                        #         continue
                        #     good = True
                        # if not good:
                        #     break
                        visited.add(tensor_hw)
                        dualps[w] = dualps.get(w, set())
                        for vlist, rc in dualpocket:
                            # if (vlist, perm_list, rc) in dualps[w]:
                            #     continue

                            dualps[w].add((vlist, rc))
                            toadd = S.One
                            # for i in range(len(vlist)):
                            #     for j in range(len(vlist[i])):
                            #         toadd *= y[i + 1] - z[vlist[i][j]]
                            #sm0 += toadd * ring(w)
                            if y is not None:
                                sm0 += vlist.polyvalue(y, z) * rc.polyvalue(y[len(vlist):], z) * ring(w)
                            else:
                                sm0 += ring(w)

        ########################333
        #check_elem = check_elem.ring.from_dict({u: check_elem[u]})
        #############3
        from sympy import expand
        diff = check_val - sm0
        # diff = DSx([]).ring.from_dict({k: sympy.sympify(vv).expand() for k, vv in diff.items() if sympy.sympify(vv).expand() != S.Zero})

        try:
            

            # assert all(expand(v99) == S.Zero for v99 in diff.values()), "Noit zor0"
            assert diff == ring.zero, "Noit zor0"
            good = True
            #assert len(visited) == len(hw_tensors), "Missed tensors"

            
            
        except AssertionError as e:
            print(f"A fail {e=} {sm0=} {u=} {v=} {w=}")
            print(f"{sm0=}")
            print(f"f{check_val=}")
            print(f"{diff=}")
            for dual in dualps.get(w, set()):
                print("=================================")
                print(dual)
                print("=================************============")
            print(f"{w}: {sm0}")
            good = False

        if good:
            print(f"Success {key} at ", time.ctime())
            with lock:
                shared_recording_dict[key] = True
        else:
            with lock:
                shared_recording_dict[key] = False
            print(f"FAIL {key} at ", time.ctime())

def is_decomposable(w):
    for i in range(1, len(w) - 1):
        coset, w_J = w.coset_decomp(*list(range(1, i + 1)), *list(range(i + 2, len(w))))
        if coset.inv == 0 and set(w_J.code[: i + 1]) != {0} and set(w_J.code[i + 2 :]) != {0}:
            return True
    return False


def recording_saver(shared_recording_dict, lock, verification_filename, stop_event, meta=None, sleep_time=10):
    """
    Background process that periodically snapshots `shared_recording_dict` and
    writes it to `verification_filename`. The lock is held only briefly to
    snapshot keys; the potentially slow file I/O and proxy lookups are done
    outside the lock so workers are not blocked.
    """
    last_verification_len_seen = -1
    while not stop_event.is_set():
        try:
            new_verification_len = len(shared_recording_dict)
        except Exception:
            new_verification_len = last_verification_len_seen

        if new_verification_len > last_verification_len_seen:
            last_verification_len_seen = new_verification_len
            print("Saving verification to", verification_filename, "with", new_verification_len, "entries at", time.ctime(), flush=True)

            # Snapshot keys quickly while holding the lock
            with lock:
                keys = list(shared_recording_dict.keys())

            # Build a stable copy outside the lock
            recording_copy = {}
            for k in keys:
                try:
                    recording_copy[k] = shared_recording_dict[k]
                except Exception:
                    # key vanished since snapshot; skip it
                    continue

            safe_save_recording(recording_copy, verification_filename, meta=meta or {})

        time.sleep(sleep_time)

    # Final save on exit: snapshot keys under lock, then copy outside
    with lock:
        keys = list(shared_recording_dict.keys())

    recording_copy = {}
    for k in keys:
        try:
            recording_copy[k] = shared_recording_dict[k]
        except Exception:
            continue

    safe_save_recording(recording_copy, verification_filename, meta=meta or {})


def queue_producer(task_queue, perms, n, num_processors, skip_id, irreducible, base_level, sep_descs, dominant_only, w0_only, shared_recording_dict, molevsagan):
    """Producer with periodic garbage collection."""
    from schubmult import Permutation, DSx, Sx
    w0 = Permutation.w0(n)
    task_count = 0
    
    for hw_tab in perms:
        if skip_id and hw_tab.inv == 0:
            continue
        if irreducible and is_decomposable(hw_tab):
            continue
        #if (not dominant_only or hw_tab.minimal_dominant_above() == hw_tab) and (not w0_only or hw_tab == w0):
        # TEMP
        if (not dominant_only or set(hw_tab.trimcode) == {1}):
            for perm in perms:
                if irreducible and is_decomposable(perm):
                    continue
                if base_level and perm[0] == 1:
                    continue
                if sep_descs:
                    if hw_tab.inv == 0 or perm.inv == 0 or max(hw_tab.descents()) <= min(perm.descents()):
                        key = (hw_tab, perm, n)
                        if shared_recording_dict.get(key) is True:
                            continue
                        task_queue.put((hw_tab, perm, n))
                else:
                    key = (hw_tab, perm, None, n)
                    if shared_recording_dict.get(key) is True:
                        continue
                    if molevsagan:
                        check_elem = DSx(hw_tab) * DSx(perm, "z")
                    else:
                        check_elem = Sx(hw_tab) * Sx(perm)
                    # Process and delete keys immediately to free memory
                    
                    task_queue.put((key, dict(check_elem)))
                    task_count += 1
                    del check_elem  # ensure cleanup
        if task_count % 100 == 0:
            gc.collect()
    
    # signal workers to exit
    for _ in range(num_processors):
        task_queue.put(None)
    

def main():
    from schubmult import DSx, Permutation, RCGraph, RootTableau, Sx, uncode  # noqa: F401

    parser = argparse.ArgumentParser(description="LR rule verification")
    parser.add_argument("n", type=int)
    parser.add_argument("num_processors", type=int)
    parser.add_argument("filename", nargs="?", help="base filename for verification output (json will be filename.json)")
    # three-state options: mutually exclusive for each bool so we can detect "not provided"
    def bool_group(pref, default=None):
        g = parser.add_mutually_exclusive_group()
        g.add_argument(f"--{pref.replace('_','-')}", dest=pref, action="store_true", help=f"enable {pref}")
        g.add_argument(f"--no-{pref.replace('_','-')}", dest=pref, action="store_false", help=f"disable {pref}")
        parser.set_defaults(**{pref: default})

    # defaults used when neither CLI nor file meta specify
    defaults = {
        "dominant_only": True,
        "w0_only": False,
        "base_level": False,
        "sep_descs": False,
        "irreducible": False,
        "skip_id": False,
        "molevsagan": False,
    }

    bool_group("dominant_only", None)
    bool_group("w0_only", None)
    bool_group("base_level", None)
    bool_group("sep_descs", None)
    bool_group("irreducible", True)
    bool_group("skip_id", None)
    bool_group("molevsagan", None)

    # Add command-line option to process subset of permutations
    parser.add_argument("--perm-range", nargs=2, type=int, help="Process perms[start:end] only")

    args = parser.parse_args()
    n = args.n
    num_processors = args.num_processors
    filename = args.filename
    verification_filename = None
    if filename:
        verification_filename = filename + ".verification"

    # Load recording and meta from JSON (if present)
    loaded_recording = {}
    loaded_meta = {}
    if verification_filename is not None:
        loaded_recording, loaded_meta = safe_load_recording(verification_filename)

    # Merge meta: CLI args override file meta; file meta overrides defaults
    meta_config = {}
    for k in defaults:
        arg_val = getattr(args, k)
        if arg_val is not None:
            meta_config[k] = bool(arg_val)
        else:
            meta_config[k] = bool(loaded_meta.get(k, defaults[k]))

    perms = Permutation.all_permutations(n)
    perms.sort(key=lambda p: (p.inv, p.trimcode))

    with Manager() as manager:
        shared_recording_dict = manager.dict()
        lock = manager.Lock()
        stop_event = Event()

        # preload records from file only (not meta)
        if verification_filename is not None and loaded_recording:
            shared_recording_dict.update(loaded_recording)

        # start saver (pass meta so saved JSON includes it)
        if verification_filename is not None:
            recording_saver_proc = Process(
                target=recording_saver, args=(shared_recording_dict, lock, verification_filename, stop_event, meta_config)
            )
            recording_saver_proc.start()

        task_queue = manager.Queue()
        # load runtime flags from meta_config
        dominant_only = meta_config["dominant_only"]
        w0_only = meta_config["w0_only"]
        base_level = meta_config["base_level"]
        sep_descs = meta_config["sep_descs"]
        irreducible = meta_config["irreducible"]
        skip_id = meta_config["skip_id"]
        molevsagan = meta_config["molevsagan"]

        

        w0 = Permutation.w0(n)
        if irreducible:
            print("Only verifying irreducible RC graphs.")
        if molevsagan:
            print("Using Molev-Sagan version of LR rule.")
        if dominant_only:
            print("Only verifying dominant multipliers.")
        if w0_only:
            print("Only verifying w0 RC.")
        if sep_descs:
            print("Testing only separated descents.")
        if skip_id:
            print("Skipping identity.")

        # queue population (unchanged logic), but skip enqueuing items already True in shared_recording_dict
        
        # Process range specified, use subset of perms
        if args.perm_range:
            start, end = args.perm_range
            perms = perms[start:end]
            print(f"Processing permutation range [{start}:{end}] ({len(perms)} perms)", flush=True)

        # Start producer process to populate queue concurrently with workers
        producer_proc = Process(
            target=queue_producer,
            args=(task_queue, perms, n, num_processors, skip_id, irreducible, base_level, sep_descs, dominant_only, w0_only, shared_recording_dict, molevsagan),
        )
        producer_proc.start()

        # Start fixed number of workers and place sentinels so they exit cleanly
        from schubmult import y, z
        workers = []
        for _ in range(num_processors):
            if molevsagan:
                p = Process(target=worker, args=(n, shared_recording_dict, lock, task_queue, y, z))
            else:
                p = Process(target=worker, args=(n, shared_recording_dict, lock, task_queue, None, None))
            p.start()
            workers.append(p)
        
        # Wait for producer to finish (it sends sentinels after all tasks)
        producer_proc.join()
        
        for p in workers:
            p.join()

        # Signal savers to exit
        stop_event.set()
        if verification_filename is not None:
            recording_saver_proc.join()

        # report results
        print_failures = True
        if print_failures:
            if any(v is False for v in shared_recording_dict.values()):
                print("Failures:")
                num_fail = 0
                for k, v in shared_recording_dict.items():
                    if v is False:
                        num_fail += 1
                        print(k)
                print(f"{num_fail} failures.")
            else:
                print("All verified successfully!")


if __name__ == "__main__":
    main()
