# LR rule verification script

from functools import cache
import gc
import json
import os
import pickle
import shutil
import sys
import time
from json import dump, load
from multiprocessing import Event, Lock, Manager, Process, cpu_count
#from schubmult.schub_lib.rc_graph_ring import tensor_to_highest_weight2
from joblib import Parallel, delayed


def reload_modules(dct, n_jobs=None):
    # from schubmult.schub_lib.rc_graph_module import RCGraph, ASx

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
        from schubmult.schub_lib.rc_graph_module import TensorModule

        with open(temp_json, "w") as f:
            json.dump(
                {repr(k): {repr(tuple(k2)): int(v2) for k2, v2 in v.value_dict.items()} if isinstance(v, TensorModule) else v for k, v in obj.items()},
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
            json.dump({repr(tuple(k)): v for k, v in obj.items()}, f)
        if os.path.exists(json_file):
            shutil.copy2(json_file, f"{json_file}.backup")
        os.replace(temp_json, json_file)
    except Exception as e:
        import traceback

        print("Error during recording JSON save:")
        traceback.print_exc()
        if os.path.exists(temp_json):
            os.remove(temp_json)


def safe_load(filename):
    pickle_file = f"{filename}.pkl"
    json_file = f"{filename}.json"
    # Try pickle first
    if os.path.exists(pickle_file):
        try:
            with open(pickle_file, "rb") as f:
                return pickle.load(f)
        except Exception as e:
            print(f"Pickle load failed: {e}")
            raise
    else:
        print(f"No pickle file {pickle_file} found.")
    # Fallback to JSON
    print("Falling back to JSON load...")
    if os.path.exists(json_file):
        try:
            with open(json_file) as f:
                loaded = json.load(f)
            print("Reloading modules from JSON backup...")
            ret = reload_modules(loaded)
            print("Done.")
            print(f"Saving as pickle to {pickle_file} for future runs...")
            safe_save(ret, filename, save_json_backup=False)
            return ret
        except Exception as e:
            print(f"JSON load failed: {e}")
            raise
    else:
        print(f"No JSON file {json_file} found.")
    return {}


def safe_load_recording(filename, Permutation):
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



# def try_lr_module(perm, length=None):
#     from schubmult.schub_lib.rc_graph_ring import RCGraphRing
#     from schubmult.schub_lib.rc_graph import RCGraph
#     from schubmult import ASx, uncode
#     ring = RCGraphRing()
#     tring = ring @ ring
#     # print(f"Starting {perm}")
#     if length is None:
#         length = len(perm.trimcode)
#     elif length < len(perm.trimcode):
#         raise ValueError("Length too short")
#     if perm.inv == 0:
#         return tring((RCGraph([()]*length),RCGraph([()]*length)))
#     if length > len(perm.trimcode):
#         mul_elem = 0
#         lower_perm = perm
#     else:
#         mul_elem = perm.trimcode[-1]
#         lower_perm = uncode(perm.trimcode[:-1])
#     elem = ASx(lower_perm, length - 1)
#     lower_module1 = try_lr_module(lower_perm, length - 1)
#     # assert isinstance(lower_module1, TensorModule), f"Not TensorModule {type(lower_module1)} {lower_perm=} {length=}"
#     #  #  # print(f"Coproducting {ASx(uncode([perm.trimcode[0]]), 1).coproduct()=}")
#     #  #  # print(ASx(uncode([perm.trimcode[0]]), 1).coproduct())
#     #  #  # print("Going for it")
#     #  #  # print(f"{type(lower_module1)=} {lower_module1=}")
#     #  #  # print(f"{type(ASx(uncode([perm.trimcode[0]]), 1).coproduct())=}")

#     cprod = tring.zero

#     for j in range(mul_elem + 1):
#         cprod += tring.ext_multiply(ring(RCGraph.one_row(j)), ring(RCGraph.one_row(mul_elem - j)))
#     #print(cprod)
#     ret_elem = lower_module1 *cprod
#     #  #  # print(f"{ret_elem=}")
#     # assert isinstance(ret_elem, TensorModule), f"Not TensorModule {type(lower_module1)} {lower_perm=} {length=}"

#     ret_elem = tensor_to_highest_weight2(tring.from_dict({k: v for k, v in ret_elem.items() if k[0].perm.bruhat_leq(perm) and k[1].perm.bruhat_leq(perm)}))

#     if length == 1:
#         return ret_elem
#     #keys = set(ret_elem.keys())
#     # print(f"{repr(keys)=} {perm=}")
#     up_elem = elem * ASx(uncode([mul_elem]),1)
#     # print(f"{up_elem=}")
#     for key, coeff in up_elem.items():
#         if key[0] != perm:
#             assert coeff == 1
#             ret_elem -= tensor_to_highest_weight2(try_lr_module(key[0], length))
#             #    ret_elem -= cff2 * tring(RCGraph.to_highest_weight_pair(rc1_bad, rc2_bad)[0])
#                 #        break
#     # print(f"Done {perm}")
#     #ret_elem = tring.from_dict({k: v for k, v in ret_elem.items() if k in keys})
#     # assert isinstance(ret_elem, TensorModule), f"Not TensorModule {type(ret_elem)} {perm.trimcode=}"
#     ret_elem = tring.from_dict({k: v for k, v in ret_elem.items() if k[0].perm.bruhat_leq(perm) and k[1].perm.bruhat_leq(perm)})
#     return ret_elem

def worker(nn, shared_recording_dict, lock, task_queue):
    from schubmult import ASx, Sx
    from sympy import pretty_print
    from schubmult import ASx, FA
    from schubmult.rings import WordBasis, SchubertBasis
    from schubmult.schub_lib.rc_graph import RCGraph
    from schubmult.schub_lib.rc_graph_ring import RCGraphRing

    ring = RCGraphRing()

    def hom(elem):
        ret = ASx.zero
        for rc, coeff in elem.items():
            ret += coeff * ASx(rc.perm, len(rc))
        return ret

    while True:
        try:
            (perm, perm2, dom_rc) = task_queue.get(timeout=2)
        except Exception:
            break  # queue empty, exit
        with lock:
            if (perm, perm2, dom_rc) in shared_recording_dict:
                if shared_recording_dict[(perm, perm2, dom_rc)] is True:
                    print(f"{(perm, perm2, dom_rc)} already verified, returning.")
                    continue
                print(f"Previous failure on {(perm, perm2, dom_rc)}, will retry.")
        result = 0
        n = len(dom_rc) + 1
        dom = RCGraph(dom_rc)
        for rc in RCGraph.all_rc_graphs(perm2, n-1):
            if rc.is_dom_perm_yamanouchi(dom.perm, perm):
                result += 1
        product = (Sx(dom.perm) * Sx(perm2))
        if result == product.get(perm, 0):
            print(f"Success {(perm, perm2, dom_rc)} at ", time.ctime())
            with lock:
                shared_recording_dict[(perm, perm2, dom_rc)] = True
        else:
            with lock:
                shared_recording_dict[(perm, perm2, dom_rc)] = False
            print(f"Warning: mismatch! {(perm, perm2, dom_rc)}: expected {product.get(perm, 0)}, got {result} at ", time.ctime())
        

def main():
    from schubmult import Permutation, uncode, Sx
    from schubmult.schub_lib.rc_graph import RCGraph

    try:
        n = int(sys.argv[1])
        filename = sys.argv[2]
        num_processors = int(sys.argv[3])
        verification_filename = filename + ".verification"
    except (IndexError, ValueError):
        print("Usage: verify_lr_rule n filename num_processors", file=sys.stderr)
        print("filename is the base file name (without extension) for saving cache results, and filename.verification is used for verification results", file=sys.stderr)
        sys.exit(1)
    cd = []
    for i in range(2 * (n - 1), 0, -2):
        cd += [i]
    
    perms = Permutation.all_permutations(n)
    perms2n = {perm for perm in Permutation.all_permutations(2 * n - 1) if perm.bruhat_leq(uncode(cd))}
    perms.sort(key=lambda p: (p.inv, p.trimcode))

    with Manager() as manager:
        shared_dict = manager.dict()
        shared_recording_dict = manager.dict()
        lock = manager.Lock()
        stop_event = Event()
        # cache_load_dict = {}
        # Load recording dict from JSON only
        loaded_recording = safe_load_recording(verification_filename, Permutation)
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
        recording_saver_proc = Process(target=recording_saver, args=(shared_recording_dict, lock, verification_filename, stop_event))
        # cache_saver_proc.start()
        recording_saver_proc.start()

        # Create task queue and fill with perms
        from schubmult.schub_lib.rc_graph import RCGraph
        task_queue = manager.Queue()
        dominant_graphs = {RCGraph.principal_rc(perm.minimal_dominant_above(), n-1) for perm in perms if perm.inv > 0 and (len(perm) - 1) <= n//2}
        for perm2 in perms:
            if perm2.inv == 0:
                continue
            for dom in dominant_graphs:
                if dom.perm.inv == 0:
                    continue
                prd = Sx(dom.perm) * Sx(perm2)
                for perm in prd.keys():
                    task_queue.put((perm, perm2, tuple(dom)))

        # Start fixed number of workers
        workers = []
        for _ in range(num_processors):
            p = Process(target=worker, args=(n, shared_recording_dict, lock, task_queue))
            p.start()
            workers.append(p)
        for p in workers:
            p.join()

        # Signal savers to exit
        stop_event.set()
        # cache_saver_proc.join()
        recording_saver_proc.join()
        print("Run finished.")
        if any(v is False for v in shared_recording_dict.values()):
            print("Failures:")
            for k, v in shared_recording_dict.items():
                if v is False:
                    print(k)
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
#     from schubmult.schub_lib.rc_graph_module import RCGraph

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
#     from schubmult.schub_lib.rc_graph_module import TensorModule

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
#     from schubmult.schub_lib.rc_graph_module import try_lr_module

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
