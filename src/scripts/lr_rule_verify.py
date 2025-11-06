# LR rule verification script

import gc
import json
import os
import pickle
import shutil
import sys
import time
from functools import cache
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

def worker(hw_tabs, nn, shared_recording_dict, lock, task_queue):
    import copy
    import sys
    from functools import cache
    from itertools import zip_longest

    import sympy
    from sympy import pretty_print

    from schubmult import CrystalGraphTensor, FreeAlgebra, Permutation, RCGraph, RCGraphRing, RootTableau, SchubertBasis, Sx


    def all_reduced_subwords(reduced_word, u):
        if u.inv > len(reduced_word):
            return set()
        if u.inv == 0:
            return {()}
        ret_set = set()
        for index in range(len(reduced_word) - 1, -1, -1):
            a = reduced_word[index]
            if a - 1 in u.descents():
                new_u = u.swap(a - 1, a)
                old_set = all_reduced_subwords(reduced_word[:index], new_u)
                for subword in old_set:
                    new_subword = (*subword, index)
                    ret_set.add(new_subword)
        return ret_set

    class MarkedInteger(int):
        pass

    hw_rc_sets = {}
    @cache
    def decompose_tensor_product(dom, u, n):
        # global hw_rc_sets
        crystals = {}
        highest_weights = set()
        perm_set = set((Sx(u)*Sx(dom.perm)).keys())
        for w in perm_set:
            if len(w) > n:
                continue
            # if not u.bruhat_leq(w):
            #     continue
            # if not dom.perm.bruhat_leq(w):
            #     continue

            # print(f"Moving on to {u=} {w=} {dom.perm=}")
            if w not in hw_rc_sets:
                hw_rc_sets[w] = set()
                for rc_w in RCGraph.all_rc_graphs(w, n - 1):
                    # pretty_print(rc_w)
                    if not rc_w.is_highest_weight:
                        continue
                    hw_rc_sets[w].add(rc_w)
            for rc_w in hw_rc_sets[w]:
                # pretty_print(rc_w)
                high_weight = rc_w.length_vector
                reduced_word = rc_w.reduced_word
                for subword in all_reduced_subwords(reduced_word, u):
                    compatible_seq = [MarkedInteger(a) if index in subword else a for index, a in enumerate(rc_w.compatible_sequence)]
                    u_tab = RootTableau.root_insert_rsk(reduced_word, compatible_seq)
                    last_inv = 1000
                    while u_tab.perm.inv < last_inv:
                        last_inv = u_tab.perm.inv
                        for box in u_tab.iter_boxes:
                            if not isinstance(u_tab[box][1], MarkedInteger):
                                u_tab_test = u_tab.delete_box(box)
                                if u_tab_test is not None:
                                    u_tab = u_tab_test
                                    break
                    if u_tab.perm.inv > u.inv:
                        # didn't make it
                        continue

                    u_tab = u_tab.rectify()
                    u_hw_rc = u_tab.rc_graph.resize(n - 1)
                    assert u_hw_rc.perm == u

                    hw_checked = set()
                    for u_tab2 in u_hw_rc.full_crystal:
                        tensor = CrystalGraphTensor(dom.rc_graph, u_tab2)
                        # print(f"{tensor=}")
                        tc_elem = tensor.to_highest_weight()[0]
                        # pretty_print(tc_elem)
                        if tc_elem in hw_checked:
                            # print("Already checked")
                            # print(f"{highest_weights=}")
                            continue
                        # needed!!!
                        if tc_elem in highest_weights:
                            # print("Already known highest weight mapped to some demazure crystal")
                            continue
                        u_tab_hw = tc_elem.factors[1]
                        # hw_checked.add(tc_elem)
                        #pretty_print(dom.rc_graph)
                        assert tc_elem.crystal_weight == tuple([a + b for a,b in zip_longest(dom.rc_graph.length_vector, u_tab_hw.length_vector, fillvalue=0)]), f"{tc_elem.crystal_weight=} vs {tuple([a + b for a,b in zip_longest(dom.rc_graph.length_vector, u_tab2.length_vector, fillvalue=0)])}"
                        high_weight_check = tuple([a for a, b in zip_longest(high_weight, tc_elem.crystal_weight, fillvalue=0)])
                        low_weight_check = tuple([a for a, b in zip_longest(rc_w.to_lowest_weight()[0].length_vector, tc_elem.crystal_weight, fillvalue=0)])
                        if tc_elem.crystal_weight == high_weight_check and tc_elem.to_lowest_weight()[0].crystal_weight == low_weight_check:
                            crystals[(rc_w, tc_elem)] = crystals.get(rc_w, 0) + 1
                            # print(f"{u=} {dom.perm=} {w=} {crystals=}")
                            highest_weights.add(tc_elem)
        return crystals

        
    ASx = FreeAlgebra(SchubertBasis)
    while True:
        try:
            (w, n) = task_queue.get(timeout=2)
        except Exception:
            break  # queue empty, exit
        with lock:
            if (w, n) in shared_recording_dict:
                if shared_recording_dict[(w, n)] is True:
                    print(f"{(w, n)} already verified, returning.")
                    continue
                print(f"Previous failure on {(w, n)}, will retry.")
        
        
        rc_w_coprods = {}
        good = False
        for hw_tab in hw_tabs:
            coprod = ASx(w, n-1).coproduct()
            for ((d, _), (u, _)) in coprod:
                if d != hw_tab.perm:
                    continue
                crystals = decompose_tensor_product(hw_tab, u, n)

                rc_ring = RCGraphRing()

                tring = rc_ring @ rc_ring

                

                for (rc_w, tc_elem), coeff in crystals.items():
                    # if not rc_w.is:
                    #     continue
                    max_len = max(len(rc_w), len(tc_elem.factors[0]), len(tc_elem.factors[1]))
                    t_elem1, t_elem2 = tc_elem.factors
                    w_rc = rc_w.resize(max_len)
                    t_elem1 = t_elem1.resize(max_len)
                    t_elem2 = t_elem2.resize(max_len)
                    if (t_elem1, t_elem2) not in rc_w_coprods.get(w_rc, tring.zero):
                        rc_w_coprods[w_rc] = rc_w_coprods.get(w_rc, tring.zero) + tring((t_elem1, t_elem2))
                    if t_elem1.perm != t_elem2.perm and (t_elem2, t_elem1) not in rc_w_coprods.get(w_rc, tring.zero):
                        rc_w_coprods[w_rc] += tring((t_elem2, t_elem1))
                    
        for rc, val in rc_w_coprods.items():
            if rc.perm != w:
                continue
            tens_ring = ASx@ASx
            check_elem = tens_ring.zero
            for (rc1, rc2), coeff in val.items():
                assert coeff == 1
                check_elem += tens_ring(((rc1.perm,len(rc1)), (rc2.perm,len(rc2))))
            diff = check_elem - coprod
            try:
                assert all(v == 0 for v in diff.values())
            except AssertionError:
                # print("A fail")
                # print(f"{diff=}")
                continue
            print(f"Coprod {rc.perm.trimcode}")
            pretty_print(rc)
            pretty_print(val)
            print("At least one success")
            good = True
            break
        #assert good, f"COMPLETE FAIL {w=}"
        if good:
            print(f"Success {(w, n)} at ", time.ctime())
            with lock:
                shared_recording_dict[(w, n)] = True
        else:
            with lock:
                shared_recording_dict[(w, n)] = False
            print(f"FAIL {(w, n)} at ", time.ctime())
    
        

def main():
    from schubmult import Permutation, RCGraph, RootTableau, Sx, uncode

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
    # perms2n = {perm for perm in Permutation.all_permutations(2 * n - 1) if perm.bruhat_leq(uncode(cd))}
    perms.sort(key=lambda p: (p.inv, p.trimcode))

    hw_tabs = set()
    for perm in perms:
        hw_tabs.update([RootTableau.from_rc_graph(rc.to_highest_weight()[0]) for rc in RCGraph.all_rc_graphs(perm, n - 1)])

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
        # dominant_graphs = {RCGraph.principal_rc(perm.minimal_dominant_above(), n-1) for perm in perms if perm.inv > 0 and (len(perm) - 1) <= n//2}
        for perm in perms:
            task_queue.put((perm, n))

        # Start fixed number of workers
        workers = []
        for _ in range(num_processors):
            p = Process(target=worker, args=(hw_tabs, n, shared_recording_dict, lock, task_queue))
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
