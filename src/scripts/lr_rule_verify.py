# LR rule verification script

import base64
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
    # print("Saving cache to", pickle_file)
    try:
        with open(temp_pickle, "wb") as f:
            pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
        if os.path.exists(pickle_file):
            shutil.copy2(pickle_file, f"{pickle_file}.backup")
        os.replace(temp_pickle, pickle_file)
    except Exception as e:
        import traceback

        # print("Error during pickle save:")
        traceback.print_exc()
        if os.path.exists(temp_pickle):
            os.remove(temp_pickle)
    # JSON backup
    if not save_json_backup:
        return
    # print("Saving JSON backup to", f"{filename}.json")
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

        # print("Error during JSON backup:")
        traceback.print_exc()
        if os.path.exists(temp_json):
            os.remove(temp_json)


def safe_save_recording(obj, filename):
    temp_json = f"{filename}.json.tmp"
    json_file = f"{filename}.json"
    try:
        # keys are permutations, values are bools
        with open(temp_json, "w") as f:
            json.dump({json_key_rep(k): v for k, v in obj.items()}, f)
        if os.path.exists(json_file):
            shutil.copy2(json_file, f"{json_file}.backup")
        os.replace(temp_json, json_file)
    except Exception as e:
        import traceback

        # print("Error during recording JSON save:")
        traceback.print_exc()
        if os.path.exists(temp_json):
            os.remove(temp_json)


def json_key_rep(k):
    """
    Serialize an arbitrary Python object `k` to a JSON-key-safe string.

    Strategy:
    - Primary: pickle the object and base64-encode the bytes, prefix with "PICKLE:".
      This preserves the original Python object on load (uses pickle.loads).
    - Fallback: if pickling fails, use "REPR:" + repr(k) and fall back to eval on load.
    """
    try:
        # use the already-imported `pickle` module from this file
        raw = pickle.dumps(k, protocol=pickle.HIGHEST_PROTOCOL)
        return "PICKLE:" + base64.b64encode(raw).decode("ascii")
    except Exception:
        # best-effort fallback
        return "REPR:" + repr(k)


def json_key_load(s):
    """
    Inverse of json_key_rep: reconstruct the original Python object from the string.

    Accepts strings produced by json_key_rep. Raises ValueError on decode/unpickle
    failure and returns the original string if no decoding rule matches.
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
            # eval is only used for the fallback REPR case
            return eval(rep)
        except Exception as e:
            raise ValueError(f"failed to eval REPR json key: {e}") from e

    # last resort: try to eval (keeps some backwards compatibility)
    try:
        return eval(s)
    except Exception:
        # give up and return the raw string
        return s


def safe_load_recording(filename):
    json_file = f"{filename}.json"
    if os.path.exists(json_file):
        try:
            with open(json_file) as f:
                loaded = json.load(f)
            # reconstruct keys as Permutation objects
            # print(f"Loaded {len(loaded)} entries from {json_file}")
            dct = {}
            for k, v in loaded.items():
                tp = json_key_load(k)
                dct[tp] = v
            return dct
        except Exception as e:
            # print(f"Recording JSON load failed: {e}")
            raise
    else:
        print(f"No recording file {json_file} found.")
    return {}


# def cache_saver(lock, filename, stop_event, sleep_time=5):
#     last_saved_results_len_seen = -1
#     while not stop_event.is_set():
#         new_len = len(shared_dict)
#         if new_len > last_saved_results_len_seen:
#             # print("Saving results to", filename, "with", new_len, "entries at ", time.ctime())
#             last_saved_results_len_seen = new_len
#             with lock:
#                 cache_copy = {k: shared_dict[k] for k in shared_dict.keys()}
#             safe_save(cache_copy, filename)
#         time.sleep(sleep_time)
#     with lock:
#         cache_copy = {k: shared_dict[k] for k in shared_dict.keys()}
#     safe_save(cache_copy, filename)
#     # print("Cache saver process exiting.")


def recording_saver(shared_recording_dict, lock, verification_filename, stop_event, sleep_time=5):
    last_verification_len_seen = -1
    while not stop_event.is_set():
        new_verification_len = len(shared_recording_dict)
        if new_verification_len > last_verification_len_seen:
            last_verification_len_seen = new_verification_len
            # print("Saving verification to ", verification_filename, " with ", new_verification_len, "entries at ", time.ctime())
            with lock:
                recording_copy = {k: shared_recording_dict[k] for k in shared_recording_dict.keys()}
            safe_save_recording(recording_copy, verification_filename)
        time.sleep(sleep_time)
    with lock:
        recording_copy = {k: shared_recording_dict[k] for k in shared_recording_dict.keys()}
    safe_save_recording(recording_copy, verification_filename)
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
#         #     # print("Saving results to", filename, "with", new_len, "entries at ", time.ctime())
#         #     last_saved_results_len_seen = new_len
#         #     with lock:
#         #         cache_copy = {k: shared_dict[k] for k in shared_dict.keys()}
#         #     safe_save(cache_copy, filename)
#         if new_verification_len > last_verification_len_seen:
#             last_verification_len_seen = new_verification_len
#             # print("Saving verification to ", verification_filename, " with ", len(shared_recording_dict), "entries at ", time.ctime())
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
#     # print("Saver process exiting.")



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

    def decompose_tensor_product(left_rc, u_rc, n):
        from schubmult import CrystalGraphTensor, RCGraphRing, RootTableau, uncode
        rc_ring = RCGraphRing()
        tring = rc_ring @ rc_ring
        # global hw_rc_sets
        crystals = {}
        #dom = RCGraph.principal_rc(uncode(left_rc.to_highest_weight()[0].length_vector), n -1)

        assert u_rc.crystal_length() == n - 1
        assert len(left_rc) == n - 1
        if u_rc.inv == 0:
            crystals[left_rc] = {CrystalGraphTensor(left_rc, u_rc)}
            return crystals
        if left_rc.inv == 0:
            crystals[u_rc] = {CrystalGraphTensor(left_rc, u_rc)}
            return crystals
        if len(u_rc) == 0:
            assert len(left_rc) == 0
            crystals[left_rc] = {CrystalGraphTensor(left_rc, u_rc)}
            return crystals
        if len(u_rc) == 1:
            assert len(left_rc) == 1
            crystals[RCGraph.one_row(len(left_rc[0]) + len(u_rc[0]))] = {CrystalGraphTensor(left_rc, u_rc)}
            return crystals
        cut_left_rc = left_rc.vertical_cut(n - 2)[0]
        cut_u = u_rc.vertical_cut(n-2)[0]

        lower_left_rc = left_rc.rowrange(1)
        lower_u = u_rc.rowrange(1)
        
        top_row_left = len(left_rc[0])
        top_row_right = len(u_rc[0])
        # print("Cutting:")
        # pretty_print(cut_left_rc)
        # pretty_print(cut_u)
        # cut_crystals = decompose_tensor_product(cut_left_rc, cut_u, n - 1)
        cut_crystals = decompose_tensor_product(lower_left_rc, lower_u, n - 1)
        # print(f"{cut_crystals=}")
        fyi = {}
        up_tensor = tring.zero
        up_rc = rc_ring.zero
        # for rc_w_cut, tensor_elems in cut_crystals.items():
        #     # up_rc =  rc_ring(rc_w_cut) * rc_ring(RCGraph.one_row(len(left_rc[-1]) + len(u_rc[-1])))
            
        #     up_rc += rc_ring(RCGraph.one_row(top_row_left + top_row_right)) * rc_ring(rc_w_cut)
        #     for t_elem in tensor_elems:
        #         # to_add =  tring(t_elem.factors) * tring((RCGraph.one_row(len(left_rc[-1])),RCGraph.one_row(len(u_rc[-1]))))
                
        #         to_add = tring((RCGraph.one_row(top_row_left),RCGraph.one_row(top_row_right))) * tring(t_elem.factors)
        #         # pretty_print(to_add)
        #         for (rc1, rc2), coeff in to_add.items():
        #             assert coeff == 1
        #             if rc1.perm != left_rc.perm or rc2.perm != u_rc.perm:
        #                 continue
        #             # tcryst = CrystalGraphTensor(rc1, rc2)
        #             # for tw in tcryst.full_crystal:
        #             #     if tw.factors[1] == u_rc:
        #             if (rc1, rc2) not in up_tensor:
        #                 up_tensor += coeff * tring((rc1, rc2))
        #             #        break
        #                 #up_tensor += coeff * tring((rc1, rc2))
            # print("up_tensor=")
            # pretty_print(up_tensor)
            # print("up_rc=")
            # pretty_print(up_rc)
        #used_tensors = set()
        for rc_w_cut, tensor_elems in cut_crystals.items():
            up_rc = rc_ring(RCGraph.one_row(top_row_left + top_row_right)) * rc_ring(rc_w_cut)
            for w_rc, coeff in up_rc.items():
                assert coeff == 1
                high_weight = w_rc.to_highest_weight()[0].crystal_weight
                low_weight = w_rc.to_lowest_weight()[0].crystal_weight
                for t_elem in tensor_elems:
                    #up_rc =  rc_ring(rc_w_cut) * rc_ring(RCGraph.one_row(len(left_rc[-1]) + len(u_rc[-1])))
                    up_tensor = tring((RCGraph.one_row(top_row_left),RCGraph.one_row(top_row_right))) * tring(t_elem.factors)
                    for (rc1, u_rc2), coeff2 in up_tensor.items():
                        assert coeff2 == 1
                        if rc1.perm != left_rc.perm or u_rc2.perm != u_rc.perm:
                            continue
                        if w_rc.perm not in (Sx(rc1.perm) * Sx(u_rc2.perm)).keys():
                            continue
                        # NOTE DOM TENSOR
                        # if rc1.perm != left_rc.perm or u_rc2.perm != u_rc.perm:
                        #     continue
                        # dom_rc = RCGraph.principal_rc(uncode(left_rc.to_highest_weight()[0].length_vector), n - 1)
                        tensor = CrystalGraphTensor(rc1, u_rc2)
                        #tensor_dom = CrystalGraphTensor(dom_rc, u_rc2)
                        tensor_hw = tensor.to_highest_weight()[0]
                        tensor_lw = tensor.to_lowest_weight()[0]
                        # if tensor_hw in used_tensors:
                        #     continue
                        #low_tensor_weight = tuple([a + b for a,b in zip(left_rc.to_lowest_weight()[0].length_vector, tensor_lw.factors[1].length_vector)])
                        low_tensor_weight = tensor_lw.crystal_weight
                        w_tab = RootTableau.from_rc_graph(w_rc)
                        u = u_rc.perm
                        if tensor_hw.crystal_weight == high_weight and low_tensor_weight == low_weight:
                            # and (u_rc2.perm.minimal_dominant_above() == u_rc2.perm or w_rc.perm.minimal_dominant_above() != w_rc.perm):
                            # for subword in all_reduced_subwords(w_rc.reduced_word, u):
                            #     # print("w_tab")
                            #     # pretty_print(w_tab)
                            #     # print(f"{w_rc.reduced_word=}")
                            #     # pretty_print(dom_tab)
                            #     roots = [w_tab.perm.right_root_at(index, word=w_rc.reduced_word) for index in subword]
                            #     grid = copy.deepcopy(w_tab._root_grid)
                            #     for box in w_tab.iter_boxes:
                            #         # print(box)
                            #         if grid[box][0] in roots:
                            #             grid[box] = (grid[box][0], MarkedInteger(grid[box][1]))
                            #     u_tab = RootTableau(grid)
                            #     last_inv = 1000
                                
                            #     while u_tab.perm.inv < last_inv:
                            #         last_inv = u_tab.perm.inv
                            #         for box in u_tab.iter_boxes_row_word_order:
                            #             if not isinstance(u_tab[box][1], MarkedInteger):
                            #                 u_tab_test = u_tab.delete_box(box)
                            #                 if u_tab_test is not None:
                            #                     u_tab = u_tab_test
                            #                     break
                            #             # else:
                            #             #     try:
                                                
                            #             #         d_tab_test = d_tab.up_jdt_slide(*box, force=True)
                            #             #         if d_tab_test is not None:
                            #             #             d_tab = d_tab_test
                            #             #     except Exception:
                            #             #         froff = False
                            #             #         # print("Couldn't up jdt")
                            #             #         # pretty_print(d_tab)
                            #             #         # print(f"{box=}")
                            #     if u_tab.perm.inv > u_rc.perm.inv:
                            #         # didn't make it
                            #         # print("No make")
                            #         # print(u_tab)
                            #         continue
                            #     u_hw_rc = u_tab.rc_graph.resize(n-1)
                            #     if u_hw_rc.perm != u:
                            #         continue
                            #     if u_hw_rc == u_rc:
                            #used_tensors.add(tensor_hw)
                            crystals[w_rc] = crystals.get(w_rc, set())
                            crystals[w_rc].add(tensor)
                    
                    # fyi[w_rc] = fyi.get(w_rc, set())
                    # fyi[w_rc].add((tensor, u_tab))
        # get_rid = 0
        # the_prins = set()
        # for ww_rc, tset in crystals.items():
        #     if ww_rc.is_principal:
        #        # get_rid.update(tset)
        #        the_prins.add(ww_rc)
        #     else:
        #        get_rid += len(tset)
        # for _ in range(get_rid):
        #     for the_prin in the_prins:
        #         try:
        #             crystals[the_prin].remove(next({c for c in crystals[the_prin] if c != CrystalGraphTensor(left_rc, u_rc)}))
        #         except Exception:
        #             break
        try:
            
            assert len(crystals) == 1
        except AssertionError:
            pass
        return crystals

        
    ASx = FreeAlgebra(SchubertBasis)
    while True:
        try:
            (u, v, n) = task_queue.get(timeout=2)
        except Exception:
            break  # queue empty, exit
        with lock:
            if (u, v, n) in shared_recording_dict:
                if shared_recording_dict[(u, v, n)] is True:
                    #print(f"{(u, v, n)} already verified, returning.")
                    continue
                print(f"Previous failure on {(u, v, n)}, will retry.")
        
        
        rc_w_coprods = {}
        good = False
        #if True:
        #hw_tab = RCGraph.principal_rc(u, n - 1).to_highest_weight()[0]
        
        mdom = Permutation.w0(n)#u.minimal_dominant_above()
        # left diff
        diff_perm = u * (~mdom)
        poly_cache = {}
        w0_prin = RCGraph.principal_rc(mdom, n)
        w0 = w0_prin.perm
        sm = Sx.zero
        rc_ring = RCGraphRing()
        for u_rc in RCGraph.all_rc_graphs(u, n):
        #if True:
            crystals0 = decompose_tensor_product(w0_prin, u_rc, n + 1)
            for rc_w_1, coeff1 in crystals0.items():
                if rc_w_1.is_principal:
                    for t_elem0 in coeff1:
                        for v_rc in RCGraph.all_rc_graphs(v, n):
                            crystals = decompose_tensor_product(w0_prin, v_rc, n + 1)
                            for rc_w, coeff in crystals.items():
                                # SANITY REMOVE PRINCIPAL
                                # rc_W principal => bijection with v_rc
                                if rc_w.is_principal:
                                    for t_elem in coeff:
                                        #sm += Sx(diff_perm * rc_w.perm)
                                        # w0 on all, same as v
                                        #if u_rc.is_principal:# and v_rc.is_principal and (w0*rc_w.perm).inv == rc_w.perm.inv - w0.inv:
                                        #if u_rc.is_principal and v_rc.is_principal:
                                            # ALWAYAS REMEMBER: rc_w principal, the crytsal is isomorphic to v_rc, with w0 weight removed
                                        for bong0 in t_elem0.full_crystal:
                                            for bong in t_elem.full_crystal:
                                                sm += Sx(bong0.factors[1].polyvalue(Sx.genset) * bong.factors[1].polyvalue(Sx.genset))

                                        # for bong0 in t_elem0.full_crystal:
                                        #     for bong in t_elem.full_crystal:
                                        #         sm += bong0.factors[1].polyvalue(Sx.genset) * bong.factors[1].polyvalue(Sx.genset)
                                    

                            
                                # take this rc_w and remove the non-u roots
                                # prin_tab = RootTableau.from_rc_graph(rc_w)
                                # dom_roots = [mdom.right_root_at(i) for i in range(mdom.inv)]
                                # u_roots = [u.right_root_at(i) for i in range(u.inv)]
                                # last_inv = 1000
                                # while prin_tab.perm.inv > u.inv + v.inv and last_inv > prin_tab.perm.inv:
                                #     last_inv = prin_tab.perm.inv
                                #     for box in prin_tab.iter_boxes:
                                #         if prin_tab[box][0] in dom_roots and prin_tab[box][0] not in u_roots:
                                #             prin_tab_test = prin_tab.delete_box(box)
                                #             if prin_tab_test is not None:
                                #                 prin_tab = prin_tab_test
                                #                 break
                                # if prin_tab.perm.inv == u.inv + v.inv:
                                #     sm += Sx(prin_tab.perm)


                                #mdom2 = uncode(
                                # for u_rc in RCGraph.all_rc_graphs(u, n - 1):
                                #     #sm += Sx(CrystalGraphTensor(u_rc, v_rc).to_lowest_weight()[0])
                                #     crystals = decompose_tensor_product(u_rc, v_rc, n)
                                #     for rc_w2, coeff2 in crystals.items():
                                #         if rc_w2.is_principal:
                                #             sm += len(coeff2) * Sx(rc_w2.perm)
                                # u_rc tensor v_rc
                                # for u_rc in RCGraph.all_rc_graphs(u, n - 1):
                                #     # decomp = decompose_tensor_product(prin2, u_rc, n)
                                #     # for rc_w2, coeff2 in decomp.items():
                                #     #     if rc_w2.is_principal:
                                #     #         minus_vec = tuple([a - b for a,b in zip(rc_w2.length_vector, prin2.length_vector)])
                                #     #         for t_elem_elem in coeff2:
                                #     #             sm += Sx(uncode(tuple([a+b for a,b in zip(minus_vec, v_rc.length_vector)])))
                                #     #             break
                                #     raise NotImplementedError("Needs work, dominant has unique subword of rc_w and this permutes and div diffs v_rc")
                    #else
                        #used_rcs.add(rc_w)


        check_elem = Sx(u) * Sx(v)
        diff = check_elem - sm
        try:
            assert diff == Sx.zero
            #print(f"Coprod {rc.perm.trimcode}")
            # pretty_print(rc)
            # pretty_print(val)
        #print("At least one success")
            good = True
        except AssertionError:
            print("A fail")
            print(f"{diff=}")
            print(f"{sm=}")
            good = False
            
        #assert good, f"COMPLETE FAIL {w=}"
        if good:
            print(f"Success {(u, v, n)} at ", time.ctime())
            with lock:
                shared_recording_dict[(u, v, n)] = True
        else:
            with lock:
                shared_recording_dict[(u, v, n)] = False
            print(f"FAIL {(u, v, n)} at ", time.ctime())
    
        

def is_decomposable(w):
    for i in range(1, len(w) - 1):
        coset, w_J = w.coset_decomp(*list(range(1, i + 1)),*list(range(i + 2, len(w))))
        if coset.inv == 0 and set(w_J.code[:i+1]) != {0} and set(w_J.code[i+2:]) != {0}:
            return True
    return False

def main():
    from schubmult import Permutation, RCGraph, RootTableau, Sx, uncode

    try:
        n = int(sys.argv[1])
        filename = sys.argv[2]
        num_processors = int(sys.argv[3])
        verification_filename = filename + ".verification"
    except (IndexError, ValueError):
        # print("Usage: verify_lr_rule n filename num_processors", file=sys.stderr)
        # print("filename is the base file name (without extension) for saving cache results, and filename.verification is used for verification results", file=sys.stderr)
        sys.exit(1)
    cd = []
    for i in range(2 * (n - 1), 0, -2):
        cd += [i]
    
    perms = Permutation.all_permutations(n)
    # perms2n = {perm for perm in Permutation.all_permutations(2 * n - 1) if perm.bruhat_leq(uncode(cd))}
    perms.sort(key=lambda p: (p.inv, p.trimcode))

    # hw_tabs = set()
    # for perm in perms:
    #     # if perm != perm.minimal_dominant_above():
    #     #     continue
    #     hw_tabs.add(RCGraph.all_rc_graphs(perm, n - 1))

    with Manager() as manager:
        shared_dict = manager.dict()
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
        #     # print(f"Loaded {len(cache_load_dict)} entries from {filename}")
        #     shared_dict.update(cache_load_dict)

        # print("Starting from ", len(shared_dict), " saved entries")
        # print("Starting from ", len(shared_recording_dict), " verified entries")
        # cache_saver_proc = Process(target=cache_saver, args=( lock, filename, stop_event))
        recording_saver_proc = Process(target=recording_saver, args=(shared_recording_dict, lock, verification_filename, stop_event))
        # cache_saver_proc.start()
        recording_saver_proc.start()

        # Create task queue and fill with perms
        from schubmult.schub_lib.rc_graph import RCGraph
        task_queue = manager.Queue()
        # dominant_graphs = {RCGraph.principal_rc(perm.minimal_dominant_above(), n-1) for perm in perms if perm.inv > 0 and (len(perm) - 1) <= n//2}
        dominant_only = False
        sep_descs = False
        indec = False
        for hw_tab in perms:
            if indec and is_decomposable(hw_tab):
                continue
            if not dominant_only or hw_tab.minimal_dominant_above() == hw_tab:
                for perm in perms:
                    if indec and is_decomposable(perm):
                        continue
                    if sep_descs:
                        if hw_tab.inv == 0 or perm.inv == 0 or max(hw_tab.descents()) <= min(perm.descents()):
                            task_queue.put((hw_tab, perm, n))
                    else:
                        task_queue.put((hw_tab, perm, n))

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
        # print("Run finished.")
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

