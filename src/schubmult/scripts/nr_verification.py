# LR rule verification script

import os
import shutil
import sys
import time
from json import dump, load
from multiprocessing import Lock, Manager, Pool, Process, cpu_count

from schubmult.rings.rc_graph_module import try_lr_module


def reload_modules(dct):
    from schubmult.rings.rc_graph_module import RCGraph

    result = {}
    for k, v in dct.items():
        key = eval(k)  # tuple
        result[key] = None
        for k2, v2 in v.items():
            g = eval(k2)
            g1, g2 = RCGraph(g[0]), RCGraph(g[1])
            if result[key] is None:
                result[key] = v2 * (g1 @ g2)
            else:
                result[key] += v2 * (g1 @ g2)
    return result


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


def saver(shared_cache_dict, shared_recording_dict, lock, max_len, filename, verification_filename, sleep_time=5):
    #last_saved_results_len_seen = 0
    last_verification_len_seen = 0
    with lock:
        #last_saved_results_len_seen = len(shared_cache_dict)
        last_verification_len_seen = len(shared_recording_dict)
    if last_verification_len_seen == max_len:
        print("Saved verification results that were loaded are complete, exiting saver at ", time.ctime())
        return
    while True:
        with lock:
            new_len = len(shared_cache_dict)
            new_verification_len = len(shared_recording_dict)
            # if new_len > last_saved_results_len_seen:
            #     print("Saving results to", filename, "with", new_len, "entries at ", time.ctime())
            #     last_saved_results_len_seen = new_len
            #     safe_save(shared_cache_dict, filename)
            if new_verification_len > last_verification_len_seen:
                last_verification_len_seen = new_verification_len
                print("Saving verification to ", verification_filename, " with ", len(shared_recording_dict), "entries at ", time.ctime())
                safe_save(shared_recording_dict, verification_filename)
            if new_verification_len >= max_len:
                print("Reached max len, exiting saver at ", time.ctime())
                return
        time.sleep(sleep_time)


def worker(args):
    shared_cache_dict, shared_recording_dict, lock, perm = args
    from schubmult import ASx
    from schubmult.rings.rc_graph_module import nonrecursive_lr_module

    with lock:
        if perm in shared_recording_dict:
            print(f"{perm} already verified, returning.")
            return  # already verified
    try_mod = nonrecursive_lr_module(perm)
    elem = try_mod.asdtype(ASx @ ASx)

    check = ASx(perm).coproduct()
    try:
        assert all(v == 0 for v in (elem - check).values())
    except AssertionError:
        print(f"Fail on {perm} at ", time.ctime())
        print(f"Module for {perm}")
        print(try_mod)
        print(f"Expected module for {perm}")
        print(try_lr_module(perm))
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
        cache_load_dict = {}
        # Load from file if it exists
        if os.path.exists(filename):
            try:
                with open(filename, "r") as f:
                    loaded = load(f)
                    if isinstance(loaded, dict):
                        cache_load_dict.update(loaded)
                        print(f"Loaded {len(loaded)} entries from {filename}")
                shared_cache_dict.update(reload_modules(cache_load_dict))
                print("Successfully reconstructed modules")
            except Exception as e:
                print(f"Could not load from {filename}: {e}")
                raise

        if os.path.exists(verification_filename):
            try:
                with open(verification_filename, "r") as f:
                    loaded = load(f)
                    if isinstance(loaded, dict):
                        shared_recording_dict.update(loaded)
                        print(f"Loaded {len(loaded)} entries from {verification_filename}")
            except Exception as e:
                print(f"Could not load from {filename}: {e}")
                raise

        print("Starting from ", len(shared_cache_dict), " saved entries")
        print("Starting from ", len(shared_recording_dict), " verified entries")
        processes = [Process(target=saver, args=(shared_cache_dict, shared_recording_dict, lock, len(perms), filename, verification_filename))]
        processes[0].start()

        # Use a process pool for workers
        pool_size = num_processors
        with Pool(processes=pool_size) as pool:
            pool.map(worker, [(shared_cache_dict, shared_recording_dict, lock, perm) for perm in perms])

        processes[0].join()

        pool.join()

        print("Verification finished.")


if __name__ == "__main__":
    main()
