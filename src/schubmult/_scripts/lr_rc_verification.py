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


def reload_modules(dct, n_jobs=None):
    # from schubmult import RCGraph, ASx

    # def reconstruct_one(k, v):
    #     key = eval(k)  # tuple
    #     result_val = None
    #     for k2, v2 in v.items():
    #         g = eval(k2)
    #         g1, g2 = RCGraph(g[0]), RCGraph(g[1])
    #         if result_val is None:
    #             result_val = v2 * (g1 @ g2)
    #         else:
    #             result_val += v2 * (g1 @ g2)
    #     return (key, result_val)

    # # Use joblib to parallelize reconstruction
    # items = list(dct.items())
    # n_jobs = n_jobs or min(8, os.cpu_count() or 1)
    # results = Parallel(n_jobs=n_jobs)(delayed(reconstruct_one)(k, v) for k, v in items)
    # return dict(results)
    return dct


def safe_save(obj, filename, save_json_backup=True):
    # Pickle save
    temp_pickle = f"{filename}.pkl.tmp"
    pickle_file = f"{filename}.pkl"
    print("Saving cache to", pickle_file)
    try:
        with open(temp_pickle, "wb") as f:
            pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
        if os.path.exists(pickle_file):
            shutil.copy2(pickle_file, f"{pickle_file}.backup")
        os.replace(temp_pickle, pickle_file)
    except Exception as e:
        import traceback

        print("Error during pickle save:")
        traceback.print_exc()
        if os.path.exists(temp_pickle):
            os.remove(temp_pickle)
    # JSON backup
    if not save_json_backup:
        return
    print("Saving JSON backup to", f"{filename}.json")
    temp_json = f"{filename}.json.tmp"
    json_file = f"{filename}.json"
    try:
        from schubmult import TensorModule

        with open(temp_json, "w") as f:
            json.dump(
                {repr(k): {repr(tuple([tuple(a) for a in k2])): int(v2) for k2, v2 in v.value_dict.items()} if isinstance(v, TensorModule) else v for k, v in obj.items()},
                f,
            )
        if os.path.exists(json_file):
            shutil.copy2(json_file, f"{json_file}.backup")
        os.replace(temp_json, json_file)
    except Exception as e:
        import traceback

        print("Error during JSON backup:")
        traceback.print_exc()
        if os.path.exists(temp_json):
            os.remove(temp_json)


def safe_save_recording(obj, filename):
    temp_json = f"{filename}.json.tmp"
    json_file = f"{filename}.json"
    try:
        # keys are permutations, values are bools
        with open(temp_json, "w") as f:
            json.dump({repr(tuple((tuple(tuple(b) for b in k[0]), tuple(k[1])))): v for k, v in obj.items()}, f)
        if os.path.exists(json_file):
            shutil.copy2(json_file, f"{json_file}.backup")
        os.replace(temp_json, json_file)
    except Exception as e:
        import traceback

        print("Error during recording JSON save:")
        traceback.print_exc()
        if os.path.exists(temp_json):
            os.remove(temp_json)



def safe_load_recording(filename):
    json_file = f"{filename}.json"
    if os.path.exists(json_file):
        try:
            with open(json_file) as f:
                loaded = json.load(f)
            # reconstruct keys as Permutation objects
            print(f"Loaded {len(loaded)} entries from {json_file}")
            dct = {}
            for k, v in loaded.items():
                tp = eval(k)
                dct[tp] = v
            return dct
        except Exception as e:
            print(f"Recording JSON load failed: {e}")
            raise
    else:
        print(f"No recording file {json_file} found.")
    return {}


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


def recording_saver(shared_recording_dict, lock, verification_filename, stop_event, sleep_time=5):
    last_verification_len_seen = -1
    while not stop_event.is_set():
        new_verification_len = len(shared_recording_dict)
        if new_verification_len > last_verification_len_seen:
            last_verification_len_seen = new_verification_len
            print("Saving verification to ", verification_filename, " with ", new_verification_len, "entries at ", time.ctime())
            with lock:
                recording_copy = {k: shared_recording_dict[k] for k in shared_recording_dict.keys()}
            safe_save_recording(recording_copy, verification_filename)
        time.sleep(sleep_time)
    with lock:
        recording_copy = {k: shared_recording_dict[k] for k in shared_recording_dict.keys()}
    safe_save_recording(recording_copy, verification_filename)
    print("Recording saver process exiting.")


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
    from schubmult import RCGraph, RCGraphRing

    rc_ring = RCGraphRing()

    while True:
        try:
            key = task_queue.get(timeout=2)
            if key is None:  # Poison pill to signal end
                break
        except Exception:
            # Timeout - check if queue is truly empty
            if task_queue.empty():
                break
            continue
            
        (g31,) = key
        rc_g = RCGraph(g31)
        g1 = rc_ring(rc_g)
        
        with lock:
            if key in shared_recording_dict:
                if shared_recording_dict[key]:
                    print(f"{key} already verified, returning.")
                    continue
                print(f"Previous failure on {g1}, will retry.")
                
        from schubmult import FreeAlgebra, SchubertBasis
        ASx = FreeAlgebra(SchubertBasis)
        gg = g1.coproduct()
        pretty_print(gg)
        tring = ASx@ASx
        g = tring.zero
        for (rc1, rc2), coeff in gg.items():
            assert coeff == 1
            g += tring(((rc1.perm, len(rc1)),(rc2.perm,len(rc2))))
        test_elem = ASx(rc_g.perm, len(rc_g)).coproduct()

        success = True
        diff = g - test_elem
        try:
            assert all(v == 0 for k, v in diff.items()), f"{tuple(diff.items())=}"
        except AssertionError as e:
            print(f"FAILURE {tuple(g1)}")
            print(f"{diff=}")
            success = False
            with lock:
                shared_recording_dict[key] = False
            continue

        with lock:
            shared_recording_dict[key] = success
        if success:
            print(f"Success {tuple(g1)} at ", time.ctime())


def main():
    from schubmult import RCGraph

    from schubmult import Permutation

    try:
        n = int(sys.argv[1])
        filename = sys.argv[2]
        num_processors = int(sys.argv[3])
        verification_filename = filename + ".verification"
    except (IndexError, ValueError):
        print("Usage: verify_lr_rule n filename num_processors", file=sys.stderr)
        print("filename is the base file name (without extension) for saving cache results, and filename.verification is used for verification results", file=sys.stderr)
        sys.exit(1)

    perms = Permutation.all_permutations(n)
    perms.sort(key=lambda p: (p.inv, p.trimcode))

    with Manager() as manager:
        shared_recording_dict = manager.dict()
        lock = manager.Lock()
        stop_event = Event()
        # cache_load_dict = {}
        # Load recording dict from JSON only
        loaded_recording = safe_load_recording(verification_filename)
        if loaded_recording:
            shared_recording_dict.update(loaded_recording)

        # Load cache dict from pickle or JSON
        # cache_load_dict = safe_load(filename)
        # if cache_load_dict:
        #     print(f"Loaded {len(cache_load_dict)} entries from {filename}")
        #     shared_dict.update(cache_load_dict)

        # print("Starting from ", len(shared_dict), " saved entries")
        print("Starting from ", len(shared_recording_dict), " verified entries")
        # cache_saver_proc = Process(target=cache_saver, args=( lock, filename, stop_event))
        # recording_saver_proc = Process(target=recording_saver, args=(shared_recording_dict, lock, verification_filename, stop_event))
        # # cache_saver_proc.start()
        # recording_saver_proc.start()

        # Create task queue and fill with perms
        task_queue = manager.Queue()
        for perm in perms:
            graphs3 = RCGraph.principal_rc(perm, n)
            # for len1 in range(len(perm.trimcode),n):
            #     for perm2 in perms:
            #         if perm2.inv == 0:
            #             continue
            #         graphs2 = RCGraph.all_rc_graphs(perm2)
            #         for len2 in range(len(perm2.trimcode),n):
            #             # for perm3 in perms:
            #             #     if perm3.inv == 0:
            #             #         continue
            #                 #graphs3 = RCGraph.all_rc_graphs(perm3)
            #                 # for len3 in range(len(perm3.trimcode), n):
            #                     for g31 in graphs3:
            #                         for g32 in graphs2:
            #                             #for g33 in graphs1:
            #                                 g1 = g31
            #                                 g2 = g32
            #                                 #g3 = g33
            task_queue.put((graphs3,))
        print(f"Enqueued {task_queue.qsize()} tasks for n={n}.")
        
        # Start fixed number of workers
        workers = []
        for _ in range(num_processors):
            p = Process(target=worker, args=(shared_recording_dict, lock, task_queue))
            p.start()
            workers.append(p)
        
        # Wait for queue to be empty
        while not task_queue.empty():
            time.sleep(1)
        
        # Send poison pills to workers
        for _ in range(num_processors):
            task_queue.put(None)
        
        # Wait for all workers to finish
        for p in workers:
            p.join()

        # Signal savers to exit
        stop_event.set()
        
        print("Run finished.")
        if any(v is False for v in shared_recording_dict.values()):
            print("Results:")
            for k, v in shared_recording_dict.items():
                if not v:
                    print(f"{k}")
        else:
            print("All verified successfully!")


if __name__ == "__main__":
    main()

# import os
# import shutil
# import sys
# import time
# from json import dump, load
# from multiprocessing import Event, Lock, Manager, Pool, Process, cpu_count

# from joblib import Parallel, delayed


# def reload_modules(dct, n_jobs=None):
#     from schubmult import RCGraph

#     def reconstruct_one(k, v):
#         key = eval(k)  # tuple
#         result_val = None
#         for k2, v2 in v.items():
#             g = eval(k2)
#             g1, g2 = RCGraph(g[0]), RCGraph(g[1])
#             if result_val is None:
#                 result_val = v2 * (g1 @ g2)
#             else:
#                 result_val += v2 * (g1 @ g2)
#         return (key, result_val)

#     # Use joblib to parallelize reconstruction
#     items = list(dct.items())
#     n_jobs = n_jobs or min(8, os.cpu_count() or 1)
#     results = Parallel(n_jobs=n_jobs)(delayed(reconstruct_one)(k, v) for k, v in items)
#     return dict(results)


# def safe_save(obj, filename):
#     from schubmult import TensorModule

#     temp_filename = f"{filename}.tmp"
#     try:
#         with open(temp_filename, "w") as f:
#             dump({repr(k): {repr(tuple([tuple(a) for a in k2])): int(v2) for k2, v2 in v.value_dict.items()} if isinstance(v, TensorModule) else v for k, v in obj.items()}, f)
#         # Atomically repla  ce the file
#         if os.path.exists(filename):
#             shutil.copy2(filename, f"{filename}.backup")  # copy, not rename
#         os.replace(temp_filename, filename)
#     except Exception as e:
#         import traceback

#         print(f"Error during safe_save:")
#         traceback.print_exc()
#         if os.path.exists(temp_filename):
#             os.remove(temp_filename)


# def saver( shared_recording_dict, lock, max_len, filename, verification_filename, stop_event, sleep_time=5):
#     last_saved_results_len_seen = 0
#     last_verification_len_seen = 0
#     with lock:
#         last_saved_results_len_seen = len(shared_dict)
#         last_verification_len_seen = len(shared_recording_dict)
#     while not stop_event.is_set():
#         cache_copy = {}
#         recording_copy = {}
#         new_len = len(shared_dict)
#         new_verification_len = len(shared_recording_dict)
#         cache_copy = {**shared_dict}
#         recording_copy = {**shared_recording_dict}
#         if new_len > last_saved_results_len_seen:
#             print("Saving results to", filename, "with", new_len, "entries at ", time.ctime())
#             last_saved_results_len_seen = new_len
#             safe_save(cache_copy, filename)
#         if new_verification_len > last_verification_len_seen:
#             last_verification_len_seen = new_verification_len
#             print("Saving verification to ", verification_filename, " with ", len(recording_copy), "entries at ", time.ctime())
#             safe_save(recording_copy, verification_filename)
#         time.sleep(sleep_time)
#     # Final save after stop
#     safe_save({**shared_dict}, filename)
#     safe_save({**shared_recording_dict}, verification_filename)
#     print("Saver process exiting.")


# def worker(args):
#      shared_recording_dict, lock, perm = args
#     from schubmult import ASx
#     from schubmult import try_lr_module

#     with lock:
#         if perm in shared_recording_dict:
#             print(f"{perm} already verified, returning.")
#             return  # already verified
#     try_mod = try_lr_module(perm, lock=lock, cache_dict=shared_dict)
#     elem = try_mod.asdtype(ASx @ ASx)

#     check = ASx(perm).coproduct()
#     try:
#         assert all(v == 0 for v in (elem - check).values())
#     except AssertionError:
#         print(f"Fail on {perm} at ", time.ctime())
#         with lock:
#             shared_recording_dict[perm] = False
#         return

#     with lock:
#         shared_recording_dict[perm] = True
#     print(f"Success {perm.trimcode} at ", time.ctime())


# def main():
#     from schubmult import Permutation

#     try:
#         n = int(sys.argv[1])
#         filename = sys.argv[2]
#         num_processors = int(sys.argv[3]) if len(sys.argv) > 3 else max(1, cpu_count() - 2)
#         verification_filename = filename + ".verification"
#     except (IndexError, ValueError):
#         print("Usage: verify_lr_rule n filename [num_processors]", file=sys.stderr)
#         print("filename is the save file for saving intermediate results, filename.verification is used for verification results", file=sys.stderr)
#         sys.exit(1)

#     perms = Permutation.all_permutations(n)
#     perms.sort(key=lambda p: (p.inv, p.trimcode))

#     with Manager() as manager:
#         shared_dict = manager.dict()
#         shared_recording_dict = manager.dict()
#         lock = manager.Lock()
#         stop_event = Event()
#         cache_load_dict = {}
#         # Load from file if it exists
#         if os.path.exists(verification_filename):
#             try:
#                 with open(verification_filename, "r") as f:
#                     loaded = load(f)
#                     if isinstance(loaded, dict):
#                         loaded = {Permutation(eval(k)): v for k, v in loaded.items()}
#                         shared_recording_dict.update(loaded)
#                         print(f"Loaded {len(loaded)} entries from {verification_filename}")
#             except Exception as e:
#                 print(f"Could not load from {filename}: {e}")
#                 raise

#         if os.path.exists(filename):
#             try:
#                 with open(filename, "r") as f:
#                     loaded = load(f)
#                     if isinstance(loaded, dict):
#                         cache_load_dict.update(loaded)
#                         print(f"Loaded {len(loaded)} entries from {filename}")
#                 print("Reconstructing modules from loaded data...")
#                 shared_dict.update(reload_modules(cache_load_dict))
#                 print("Successfully reconstructed modules")
#             except Exception as e:
#                 print(f"Could not load from {filename}: {e}")
#                 raise

#         print("Starting from ", len(shared_dict), " saved entries")
#         print("Starting from ", len(shared_recording_dict), " verified entries")
#         saver_proc = Process(target=saver, args=( shared_recording_dict, lock, len(perms), filename, verification_filename, stop_event))
#         saver_proc.start()

#         # Use a process pool for workers
#         pool_size = num_processors
#         with Pool(processes=pool_size) as pool:
#             pool.map(worker, [( shared_recording_dict, lock, perm) for perm in perms])

#         # Signal saver to exit
#         stop_event.set()
#         saver_proc.join()
#         print("Verification finished.")


# if __name__ == "__main__":
#     main()
