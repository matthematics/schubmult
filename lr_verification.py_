# LR rule verification script

import gc
import json
import os
import pickle
import shutil
import sys
import time
from json import dump, load
from math import perm
from multiprocessing import Event, Lock, Manager, Process, cpu_count

from joblib import Parallel, delayed
from sympy import pretty_print

from schubmult import CrystalGraphTensor
import argparse
import base64


def safe_save_recording(obj, filename, meta=None):
    temp_json = f"{filename}.json.tmp"
    json_file = f"{filename}.json"
    try:
        payload = {
            "meta": meta or {},
            "records": {json_key_rep(k): v for k, v in obj.items()},
        }
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
    Load a verification JSON file.

    Backwards compatible with the older format (flat mapping): returns (records, meta).
    New format stores {"meta": {...}, "records": {...}}.
    """
    json_file = f"{filename}.json"
    if os.path.exists(json_file):
        try:
            with open(json_file) as f:
                loaded = json.load(f)
            # new-style wrapper
            if isinstance(loaded, dict) and "records" in loaded:
                raw_records = loaded.get("records", {})
                meta = loaded.get("meta", {}) or {}
            else:
                # legacy: treat whole file as records mapping
                raw_records = loaded
                meta = {}
            dct = {}
            for k, v in raw_records.items():
                tp = json_key_load(k)
                dct[tp] = v
            return dct, meta
        except Exception:
            raise
    else:
        print(f"No recording file {json_file} found.")
    return {}, {}


# def cache_saver(lock, filename, stop_event, sleep_time=5):
#     last_saved_results_len_seen = -1
#     while not stop_event.is_set():
#         new_len = len(shared_dict)
#         if new_len > last_saved_results_len_seen:
#             print("Saving results to", filename, "with", new_len, "entries at ", time.ctime())
#             last_saved_results_len_seen = new_len
#             with lock:
#                 cache_copy = {k: shared_dict[k] for k in shared_dict.keys()}
#             safe_save(cache_copy, filename)
#         time.sleep(sleep_time)
#     with lock:
#         cache_copy = {k: shared_dict[k] for k in shared_dict.keys()}
#     safe_save(cache_copy, filename)
#     print("Cache saver process exiting.")


def recording_saver(shared_recording_dict, lock, verification_filename, stop_event, meta=None, sleep_time=10):
    last_verification_len_seen = -1
    while not stop_event.is_set():
        new_verification_len = len(shared_recording_dict)
        if new_verification_len > last_verification_len_seen:
            last_verification_len_seen = new_verification_len
            print("Saving verification to ", verification_filename, " with ", new_verification_len, "entries at ", time.ctime())

            # Snapshot keys quickly while holding the lock, then build the copy outside.
            with lock:
                keys = list(shared_recording_dict.keys())

            recording_copy = {}
            for k in keys:
                try:
                    recording_copy[k] = shared_recording_dict[k]
                except Exception:
                    # key vanished since snapshot; skip it
                    continue

            safe_save_recording(recording_copy, verification_filename, meta=meta or {})
        time.sleep(sleep_time)

    # final save on exit: snapshot keys under lock, then copy outside
    with lock:
        keys = list(shared_recording_dict.keys())

    recording_copy = {}
    for k in keys:
        try:
            recording_copy[k] = shared_recording_dict[k]
        except Exception:
            continue

    safe_save_recording(recording_copy, verification_filename, meta=meta or {})
    # print("Recording saver process exiting.")


# def saver(shared_recording_dict, lock, filename, verification_filename, stop_event, sleep_time=5):
#     # last_saved_results_len_seen = -1
#     last_verification_len_seen = -1
#     new_len = 0
#     new_verification_len = 0
#     while not stop_event.is_set():
#         # print("Im in ur saver")
#         # new_len = len(shared_dict)
#         new_verification_len = len(shared_recording_dict)
#         # if new_len > last_saved_results_len_seen:
#         #     print("Saving results to", filename, "with", new_len, "entries at ", time.ctime())
#         #     last_saved_results_len_seen = new_len
#         #     with lock:
#         #         cache_copy = {k: shared_dict[k] for k in shared_dict.keys()}
#         #     safe_save(cache_copy, filename)
#         if new_verification_len > last_verification_len_seen:
#             last_verification_len_seen = new_verification_len
#             print("Saving verification to ", verification_filename, " with ", len(shared_recording_dict), "entries at ", time.ctime())
#             with lock:
#                 recording_copy = {k: shared_recording_dict[k] for k in shared_recording_dict.keys()}
#             safe_save_recording(recording_copy, verification_filename)
#         # print("me sleep")
#         time.sleep(sleep_time)
#     with lock:
#         cache_copy = {k: shared_dict[k] for k in shared_dict.keys()}
#         recording_copy = {k: shared_recording_dict[k] for k in shared_recording_dict.keys()}
#     safe_save(cache_copy, filename)
#     safe_save_recording(recording_copy, verification_filename)
#     print("Saver process exiting.")


def worker(shared_recording_dict, lock, task_queue):
    from schubmult import RCGraphRing, Sx

    rc_ring = RCGraphRing()
    while True:
        try:
            key = task_queue.get(timeout=2)
            # (g31, g32, (len1, len2)) = key
            g11, g22, g33, (len1, len2, len3) = key
            g1 = rc_ring(g11.resize(len1))
            g2 = rc_ring(g22.resize(len2))
            g3 = rc_ring(g33.resize(len3))

        except Exception:
            break  # queue empty, exit
        with lock:
            if key in shared_recording_dict:
                if shared_recording_dict[key]:
                    print(f"{key} already verified, returning.")
                    continue
                print(f"Previous failure on {(g1, g2)}, will retry.")
        g = g1 * (g2 * g3)
        g_ = (g1 * g2) * g3
        diff = g - g_
        try:
            assert all(v == 0 for k, v in diff.items()), f"{tuple(diff.items())=}"
            success = True
        except AssertionError as e:
            print("FAILURE")
            print(e)
            print(f"{g=}")
            print(f"{g_=}")
            success = False

        with lock:
            shared_recording_dict[key] = success
        if success:
            print(f"Success {key} at ", time.ctime())

        del g1
        del g2
        del g3
        del g
        del g_
        del diff
        gc.collect()


def json_key_rep(k):
    """
    Serialize an arbitrary Python object `k` to a JSON-key-safe string.
    Uses pickle + base64 with a 'PICKLE:' prefix.
    """
    try:
        raw = pickle.dumps(k, protocol=pickle.HIGHEST_PROTOCOL)
        return "PICKLE:" + base64.b64encode(raw).decode("ascii")
    except Exception:
        # Fallback: use repr (best-effort; json_key_load will try eval)
        return "REPR:" + repr(k)


def json_key_load(s):
    """
    Inverse of json_key_rep: reconstruct the original Python object.
    Accepts strings produced by json_key_rep. Raises ValueError on failure.
    """
    if not isinstance(s, str):
        raise TypeError("json_key_load expects a string")

    if s.startswith("PICKLE:"):
        payload = s[len("PICKLE:") :]
        try:
            raw = base64.b64decode(payload.encode("ascii"))
            return pickle.loads(raw)
        except Exception as e:
            raise ValueError(f"failed to unpickle json key: {e}") from e

    if s.startswith("REPR:"):
        rep = s[len("REPR:") :]
        try:
            return eval(rep)
        except Exception as e:
            raise ValueError(f"failed to eval REPR json key: {e}") from e

    # last resort: try eval, then return raw string
    try:
        return eval(s)
    except Exception:
        return s


def main():
    from schubmult import DSx, Permutation, RCGraph, RootTableau, Sx, uncode  # noqa: F401

    parser = argparse.ArgumentParser(description="LR rule verification")
    parser.add_argument("n", type=int)
    parser.add_argument("num_processors", type=int)
    parser.add_argument("filename", nargs="?", help="base filename for verification output (json will be filename.json)")
    # three-state options: mutually exclusive pair for each boolean so we can detect "not provided"
    def bool_group(pref):
        g = parser.add_mutually_exclusive_group()
        g.add_argument(f"--{pref.replace('_','-')}", dest=pref, action="store_true", help=f"enable {pref}")
        g.add_argument(f"--no-{pref.replace('_','-')}", dest=pref, action="store_false", help=f"disable {pref}")
        parser.set_defaults(**{pref: None})

    # define the flags we want to control
    parser.set_defaults()  # ensure parser has defaults dict
    bool_group("dominant_only")
    bool_group("w0_only")
    bool_group("base_level")
    bool_group("sep_descs")
    bool_group("irreducible")
    bool_group("skip_id")

    args = parser.parse_args()
    n = args.n
    num_processors = args.num_processors
    filename = args.filename
    verification_filename = None
    if filename:
        verification_filename = filename + ".verification"

    # default behavior if not specified anywhere
    defaults = {
        "dominant_only": True,
        "w0_only": False,
        "base_level": False,
        "sep_descs": False,
        "irreducible": True,
        "skip_id": True,
    }

    # Load recording dict and meta from JSON only
    loaded_recording = {}
    loaded_meta = {}
    if verification_filename is not None:
        loaded_recording, loaded_meta = safe_load_recording(verification_filename)
        if loaded_recording:
            # merge loaded records into shared dict below
            pass

    # merge meta: CLI args override file meta; file meta overrides defaults
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

        # Load recording dict from JSON only (records)
        if verification_filename is not None and loaded_recording:
            shared_recording_dict.update(loaded_recording)

            # pass meta_config to saver process so JSON meta reflects current config
            recording_saver_proc = Process(
                target=recording_saver, args=(shared_recording_dict, lock, verification_filename, stop_event, meta_config)
            )
            recording_saver_proc.start()
        elif verification_filename is not None:
            # start saver even if no initial records so it will save when work produces entries
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

        w0 = Permutation.w0(n)
        if irreducible:
            print("Only verifying irreducible RC graphs.")
        if dominant_only:
            print("Only verifying dominant multipliers.")
        if w0_only:
            print("Only verifying w0 RC.")
        if sep_descs:
            print("Testing only separated descents.")
        if skip_id:
            print("Skipping identity.")


        for hw_tab in perms:
            if skip_id and hw_tab.inv == 0:
                continue
            if irreducible and is_decomposable(hw_tab):
                continue
            if (not dominant_only or hw_tab.minimal_dominant_above() == hw_tab) and (not w0_only or hw_tab == w0):
                for perm in perms:
                    if irreducible and is_decomposable(perm):
                        continue
                    if base_level and perm[0] == 1:
                        continue
                    if sep_descs:
                        if hw_tab.inv == 0 or perm.inv == 0 or max(hw_tab.descents()) <= min(perm.descents()):
                            # when using the sep_descs path the enqueued tuple shape is different;
                            # check existing recording for that exact key shape before enqueueing.
                            key = (hw_tab, perm, n)
                            if shared_recording_dict.get(key) is True:
                                continue
                            task_queue.put((hw_tab, perm, n))
                    else:
                        check_elem = DSx(hw_tab) * DSx(perm, "z")
                        keys = set(check_elem.keys())
                        for w in keys:
                            val = check_elem[w]
                            key = (hw_tab, perm, w, n)
                            # Skip queuing if this key was already verified True in the loaded recording.
                            if shared_recording_dict.get(key) is True:
                                # already verified successfully — skip
                                continue
                            task_queue.put((key, val))
                        for w in keys:
                            del check_elem[w]
                        #gc.collect()

        # Start fixed number of workers
        workers = []

        for _ in range(num_processors):
            task_queue.put(None  # sentinel to signal worker exit
            p = Process(target=worker, args=(n, shared_recording_dict, lock, task_queue))
            p.start()
            workers.append(p)
        for p in workers:
            p.join()

        # Signal savers to exit
        stop_event.set()
        if verification_filename is not None:
            recording_saver_proc.join()
        # print("Run finished.")
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

