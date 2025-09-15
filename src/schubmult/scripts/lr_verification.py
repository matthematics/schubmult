# LR rule verification script

import os
import sys
import time
from multiprocessing import Lock, Manager, Pool, Process, cpu_count
from pickle import dump, load


def saver(shared_cache_dict, shared_recording_dict, lock, max_len, filename, verification_filename, sleep_time=5):
    last_saved_results_len_seen = 0
    last_verification_len_seen = 0
    with lock:
        last_saved_results_len_seen = len(shared_cache_dict)
        last_verification_len_seen = len(shared_recording_dict)
    if last_verification_len_seen == max_len:
        print("Saved verification results that were loaded are complete, exiting saver at ", time.ctime())
        return
    while True:
        with lock:
            new_len = len(shared_cache_dict)
            new_verification_len = len(shared_recording_dict)
            if new_len > last_saved_results_len_seen:
                print("Saving results to", filename, "with", new_len, "entries at ", time.ctime())
                last_saved_results_len_seen = new_len
                with open(filename, "wb") as f:
                    dump(dict(shared_cache_dict), f)
            if new_verification_len > last_verification_len_seen:
                last_verification_len_seen = new_verification_len
                print("Saving verification to ", verification_filename, " with ", len(shared_recording_dict), "entries at ", time.ctime())
                with open(verification_filename, "wb") as f:
                    dump(dict(shared_recording_dict), f)
            if new_verification_len >= max_len:
                print("Reached max len, exiting saver at ", time.ctime())
                return
        time.sleep(sleep_time)


def worker(args):
    shared_cache_dict, shared_recording_dict, lock, perm = args
    from schubmult import ASx
    from schubmult.rings.rc_graph_module import try_lr_module

    try_mod = None
    with lock:
        if perm in shared_recording_dict:
            print(f"{perm} already verified, returning.")
            return  # already verified
        if perm in shared_cache_dict:
            try_mod = shared_cache_dict.get(perm, None)
    if try_mod is None:
        try_mod = try_lr_module(perm)
        with lock:
            shared_cache_dict[perm] = try_mod
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

    n = int(sys.argv[1])
    filename = sys.argv[2]
    num_processors = int(sys.argv[3]) if len(sys.argv) > 3 else max(1, cpu_count() - 2)
    verification_filename = filename + ".verification"

    perms = Permutation.all_permutations(n)
    perms.sort(key=lambda p: (p.inv, p.trimcode))

    with Manager() as manager:
        shared_cache_dict = manager.dict()
        shared_recording_dict = manager.dict()
        lock = manager.Lock()

        # Load from file if it exists
        if os.path.exists(filename):
            try:
                with open(filename, "rb") as f:
                    loaded = load(f)
                    if isinstance(loaded, dict):
                        shared_cache_dict.update(loaded)
                        print(f"Loaded {len(loaded)} entries from {filename}")
            except Exception as e:
                print(f"Could not load from {filename}: {e}")

        if os.path.exists(verification_filename):
            try:
                with open(verification_filename, "rb") as f:
                    loaded = load(f)
                    if isinstance(loaded, dict):
                        shared_recording_dict.update(loaded)
                        print(f"Loaded {len(loaded)} entries from {filename}")
            except Exception as e:
                print(f"Could not load from {filename}: {e}")

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
