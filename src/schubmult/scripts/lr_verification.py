# LR rule verification script

import os
import shutil
import sys
import time
from json import dump, load
from multiprocessing import Event, Lock, Manager, Pool, Process, cpu_count

from joblib import Parallel, delayed


def reload_modules(dct, n_jobs=None):
    from schubmult.rings.rc_graph_module import RCGraph

    def reconstruct_one(k, v):
        key = eval(k)  # tuple
        result_val = None
        for k2, v2 in v.items():
            g = eval(k2)
            g1, g2 = RCGraph(g[0]), RCGraph(g[1])
            if result_val is None:
                result_val = v2 * (g1 @ g2)
            else:
                result_val += v2 * (g1 @ g2)
        return (key, result_val)

    # Use joblib to parallelize reconstruction
    items = list(dct.items())
    n_jobs = n_jobs or min(8, os.cpu_count() or 1)
    results = Parallel(n_jobs=n_jobs)(delayed(reconstruct_one)(k, v) for k, v in items)
    return dict(results)


def safe_save(obj, filename):
    from schubmult.rings.rc_graph_module import TensorModule

    temp_filename = f"{filename}.tmp"
    try:
        with open(temp_filename, "w") as f:
            dump({repr(k): {repr(tuple([tuple(a) for a in k2])): int(v2) for k2, v2 in v.value_dict.items()} if isinstance(v, TensorModule) else v for k, v in obj.items()}, f)
        # Atomically repla  ce the file
        if os.path.exists(filename):
            shutil.copy2(filename, f"{filename}.backup")  # copy, not rename
        os.replace(temp_filename, filename)
    except Exception as e:
        import traceback

        print(f"Error during safe_save:")
        traceback.print_exc()
        if os.path.exists(temp_filename):
            os.remove(temp_filename)


def saver(shared_cache_dict, shared_recording_dict, lock, max_len, filename, verification_filename, stop_event, sleep_time=5):
    last_saved_results_len_seen = 0
    last_verification_len_seen = 0
    with lock:
        last_saved_results_len_seen = len(shared_cache_dict)
        last_verification_len_seen = len(shared_recording_dict)
    while not stop_event.is_set():
        cache_copy = {}
        recording_copy = {}
        new_len = len(shared_cache_dict)
        new_verification_len = len(shared_recording_dict)
        cache_copy = {**shared_cache_dict}
        recording_copy = {**shared_recording_dict}
        if new_len > last_saved_results_len_seen:
            print("Saving results to", filename, "with", new_len, "entries at ", time.ctime())
            last_saved_results_len_seen = new_len
            safe_save(cache_copy, filename)
        if new_verification_len > last_verification_len_seen:
            last_verification_len_seen = new_verification_len
            print("Saving verification to ", verification_filename, " with ", len(recording_copy), "entries at ", time.ctime())
            safe_save(recording_copy, verification_filename)
        time.sleep(sleep_time)
    # Final save after stop
    safe_save({**shared_cache_dict}, filename)
    safe_save({**shared_recording_dict}, verification_filename)
    print("Saver process exiting.")


def worker(args):
    shared_cache_dict, shared_recording_dict, lock, perm = args
    from schubmult import ASx
    from schubmult.rings.rc_graph_module import try_lr_module_cache

    with lock:
        if perm in shared_recording_dict:
            print(f"{perm} already verified, returning.")
            return  # already verified
    try_mod = try_lr_module_cache(perm, lock=lock, cache_dict=shared_cache_dict)
    elem = try_mod.asdtype(ASx @ ASx)

    check = ASx(perm).coproduct()
    try:
        assert all(v == 0 for v in (elem - check).values())
    except AssertionError:
        print(f"Fail on {perm} at ", time.ctime())
        with lock:
            shared_recording_dict[perm] = False
        return

    with lock:
        shared_recording_dict[perm] = True
    print(f"Success {perm.trimcode} at ", time.ctime())


def main():
    from schubmult import Permutation

    try:
        n = int(sys.argv[1])
        filename = sys.argv[2]
        num_processors = int(sys.argv[3]) if len(sys.argv) > 3 else max(1, cpu_count() - 2)
        verification_filename = filename + ".verification"
    except (IndexError, ValueError):
        print("Usage: verify_lr_rule n filename [num_processors]", file=sys.stderr)
        print("filename is the save file for saving intermediate results, filename.verification is used for verification results", file=sys.stderr)
        sys.exit(1)

    perms = Permutation.all_permutations(n)
    perms.sort(key=lambda p: (p.inv, p.trimcode))

    with Manager() as manager:
        shared_cache_dict = manager.dict()
        shared_recording_dict = manager.dict()
        lock = manager.Lock()
        stop_event = manager.Event()
        cache_load_dict = {}
        # Load from file if it exists
        if os.path.exists(verification_filename):
            try:
                with open(verification_filename, "r") as f:
                    loaded = load(f)
                    if isinstance(loaded, dict):
                        # assert type(eval(next(iter(loaded.keys())))) == tuple, f"{type(eval(next(iter(loaded.keys()))))=}"
                        loaded = {Permutation(eval(k)): v for k, v in loaded.items()}
                        shared_recording_dict.update(loaded)
                        print(f"Loaded {len(loaded)} entries from {verification_filename}")
            except Exception as e:
                print(f"Could not load from {filename}: {e}")
                raise

        if os.path.exists(filename):
            try:
                with open(filename, "r") as f:
                    loaded = load(f)
                    if isinstance(loaded, dict):
                        cache_load_dict.update(loaded)
                        print(f"Loaded {len(loaded)} entries from {filename}")
                print("Reconstructing modules from loaded data...")
                shared_cache_dict.update(reload_modules(cache_load_dict))
                print("Successfully reconstructed modules")
            except Exception as e:
                print(f"Could not load from {filename}: {e}")
                raise

        print("Starting from ", len(shared_cache_dict), " saved entries")
        print("Starting from ", len(shared_recording_dict), " verified entries")
        saver_proc = Process(target=saver, args=(shared_cache_dict, shared_recording_dict, lock, len(perms), filename, verification_filename, stop_event))
        saver_proc.start()

        # Use a process pool for workers
        pool_size = num_processors
        with Pool(processes=pool_size) as pool:
            pool.map(worker, [(shared_cache_dict, shared_recording_dict, lock, perm) for perm in perms])

        # Signal saver to exit
        stop_event.set()
        saver_proc.join()
        print("Verification finished.")


if __name__ == "__main__":
    main()
