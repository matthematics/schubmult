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

def divdiff_rc_desc(the_rc, desc):
    from sympy import pretty_print

    from schubmult import RCGraph
    ret = set()
    rc, row = the_rc.exchange_property(desc, return_row=True)
    rc = rc.normalize()
    if row != desc:
        return ret
    if rc.raising_operator(desc) is not None:
        return ret
    ret.add(rc)
    while rc is not None:
        rc = rc.lowering_operator(desc)
        if rc is not None:
            ret.add(rc)
                # if the_rc0.is_lowest_weight:
                #     ret.update(addup.crystal_beneath)
    return ret

def divdiff_tensor_rc_desc(the_rc1, the_rc2, desc):
    from sympy import pretty_print

    from schubmult import CrystalGraphTensor, RCGraph
    ret = set()

    try:
        rc, row = the_rc1.exchange_property(desc, return_row=True)
        rc = rc.normalize()
        
        rc2 = the_rc2#.weight_reflection(desc)
        max_size = max(len(rc), len(rc2))
        tens = CrystalGraphTensor(rc.resize(max_size), rc2.resize(max_size))
        if tens.raising_operator(desc) is None:        
            ret.add(tens)
            while tens is not None:
                tens = tens.lowering_operator(desc)
                if tens is not None:
                    ret.add(tens)
    except ValueError:
        pass

    try:
        rc, row = the_rc2.exchange_property(desc, return_row=True)
        rc = rc.normalize()
        max_size = max(len(the_rc1), len(rc))
        tens = CrystalGraphTensor(the_rc1.resize(max_size), rc.resize(max_size))
        if tens.raising_operator(desc) is None:
            ret.add(tens)
            while tens is not None:
                tens = tens.lowering_operator(desc)
                if tens is not None:
                    ret.add(tens)
    except ValueError:
        pass

                # if the_rc0.is_lowest_weight:
                #     ret.update(addup.crystal_beneath)
    return ret


def divdiffable_rc(v_rc, u):
    from schubmult import RCGraph
    # if not v_rc.is_highest_weight:
    #     raise ValueError("Crystal")
    #     vr_hw, raise_seq = v_rc.to_highest_weight()
    #     hw_set = divdiffable_rc(vr_hw, u)
    #     ret = set()
    #     for vv_hw in hw_set:
    #         ret.add(vv_hw.reverse_raise_seq(raise_seq))
    v = v_rc.perm
    perm2 = v * (~u)
    if perm2.inv != v.inv - u.inv:
        return set()
    # return perm2
    ret = {v_rc}
    working_perm = u
    while working_perm.inv > 0:
        working_set = set()
        desc = max(working_perm.descents()) + 1
        working_perm = working_perm.swap(desc - 1, desc)
        for the_rc in ret:
            assert desc-1 in the_rc.perm.descents()
            working_set.update(divdiff_rc_desc(the_rc, desc))

        ret = working_set
    return ret


# def dualpieri(mu, v, w):
#     # logger.debug(f"dualpieri {mu=} {v=} {w=}")
#     lm = code(~mu)
#     cn1w = code(~w)
#     while len(lm) > 0 and lm[-1] == 0:
#         lm.pop()
#     while len(cn1w) > 0 and cn1w[-1] == 0:
#         cn1w.pop()
#     if len(cn1w) < len(lm):
#         return []
#     for i in range(len(lm)):
#         if lm[i] > cn1w[i]:
#             return []
#     c = Permutation([1, 2])
#     for i in range(len(lm), len(cn1w)):
#         c = cycle(i - len(lm) + 1, cn1w[i]) * c
#     # c = permtrim(c)
#     # logger.debug("Recording line number")
#     res = [[[], v]]
#     # logger.debug(f"{v=} {type(v)=}")
#     for i in range(len(lm)):
#         # logger.debug(f"{res=}")
#         res2 = []
#         for vlist, vplist in res:
#             vp = vplist
#             vpl = divdiffable(vp, cycle(lm[i] + 1, cn1w[i] - lm[i]))
#             # logger.debug(f"{vpl=} {type(vpl)=}")
#             if len(vpl) == 0:
#                 continue
#             vl = pull_out_var(lm[i] + 1, vpl)
#             # logger.debug(f"{vl=}")
#             for pw, vpl2 in vl:
#                 res2 += [[[*vlist, pw], vpl2]]
#         res = res2
#     if len(lm) == len(cn1w):
#         return res
#     res2 = []
#     for vlist, vplist in res:
#         vp = vplist
#         vpl = divdiffable(vp, c)
#         if len(vpl) == 0:
#             continue
#         res2 += [[vlist, vpl]]
#     # logger.debug(f"{res2=}")
#     return res2

def divdiff_rc_ring(rc_ring_elem, u):
    result = rc_ring_elem.ring.zero
    for rc, coeff in rc_ring_elem.items():
        rc_set = divdiffable_rc(rc, u)
        for rc0 in rc_set:
            result += coeff * rc_ring_elem.ring(rc0)
    return result

def divdiff_rc_tensor_ring(rc_ring_elem, desc):
    result = rc_ring_elem.ring.zero
    for (rc1, rc2), coeff in rc_ring_elem.items():
        rc_set = divdiff_tensor_rc_desc(rc1, rc2, desc)
        for rc0 in rc_set:
            result += coeff * rc_ring_elem.ring(rc0.factors)
    return result

def monk_mul(rc_elem, k):
    from schubmult import Permutation, RCGraph
    if k == 0:
        return rc_elem
    result = rc_elem.ring.zero
    for rc, coeff in rc_elem.items():
        index = k
        while index > 0:
            if index == k:
                start = rc.kogan_kumar_insert(k, [index])
            elif start.raising_operator(index) is not None:
                start = start.raising_operator(index)
            else:
                start = rc.kogan_kumar_insert(k, [index])
            result += coeff * rc_elem.ring(start)
            index -= 1
    return result





# def h_mul(rc_elem, p, k):
#     for rc
#     return result

# def monk_formula(rc, k):
#     from schubmult import RCGraphRing
#     rc_ring = RCGraphRing()
#     if k == 0:
#         return rc_ring(rc)
#     if k < 0:
#         return rc_ring.zero
#     if k == 1:
#         rc2 = rc.kogan_kumar_insert([1], 1)
#         return rc_ring(rc2)
#     if k == 2:
#         if rc.perm[0] < rc.perm[1]:
#             return divdiff_rc_ring(rc_ring(rc.kogan_kumar_insert([1, 1], 1)),Permutation.ref_product(1))
#         return monk_formula()
    


def dualpieri(mu, v_rc, w):
    from sympy import pretty_print

    from schubmult import Permutation, RCGraph, pull_out_var
    if mu.inv == 0:
        #ret = set()
        return set({((), rc) for rc in divdiffable_rc(v_rc, w)})
        
    cycle = Permutation.cycle
    lm = (~mu).trimcode
    cn1w = (~w).trimcode
    if len(cn1w) < len(lm):
        return set()
    for i in range(len(lm)):
        if lm[i] > cn1w[i]:
            return set()
    #lm = [*lm, *([0] * (len(cn1w) - len(lm)))]
    c = Permutation([])
    for i in range(len(lm), len(cn1w)):
        c = cycle(i - len(lm) + 1, cn1w[i]) * c 
    # c = permtrim(c)
    # logger.debug("Recording line number")
    res = {((), v_rc)}
    # logger.debug(f"{v=} {type(v)=}")
    for i in range(len(lm)):
        # logger.debug(f"{res=}")
        res2 = set()
        for vlist, v_rc_0 in res:
            vp = v_rc_0
            
            vpl_list = divdiffable_rc(vp, cycle(lm[i] + 1, cn1w[i] - lm[i]))
            # logger.debug(f"{vpl=} {type(vpl)=}")
            if len(vpl_list) == 0:
                continue
            for vpl in vpl_list:
                vl = pull_out_var(lm[i] + 1, vpl.perm)
                if lm[i] >= len(vpl):
                    res2.add(((*vlist, ()), vpl))
                    continue
                # print(f"Begin rc pulling out lm[i] + 1={lm[i] + 1}")
                # pretty_print(vp)
                # print("Before")
                # pretty_print(vpl_new)
                # for ref_spot_i in range(lm[i] + 1, len(vpl) - 1):
                #     # print(f"Before ref {ref_spot_i=}")
                #     # pretty_print(vpl_new)
                #     vpl_new = vpl_new.weight_reflection(ref_spot_i)
                #     # print(f"Reflected {ref_spot_i=}")
                #     # pretty_print(vpl_new)
                #     # if vpl_new.perm != working_perm:
                #     #     raise ValueError
                # vpl_new = vpl_new.vertical_cut(len(vpl.perm.trimcode))[0].resize(len(vpl))
                # vpl_new = vpl_new.vertical_cut(len(vpl.perm.trimcode) - 1)[0] if lm[i] + 1 <= len(vpl.perm.trimcode) else vpl_new
                rc_to_match = vpl.vertical_cut(lm[i])[0].to_highest_weight()[0]
                weight = tuple([*vpl.length_vector[:lm[i]],*vpl.length_vector[lm[i]+1:]])
                for pw, vpl2 in vl:
                    if len(vpl[lm[i]]) != len(pw):
                        continue
                    for vpl2_rc in RCGraph.all_rc_graphs(vpl2, len(vpl) - 1, weight=weight):
                        if vpl2_rc.rowrange(lm[i] + 1) == vpl.rowrange(lm[i] + 1):
                            if vpl2_rc.vertical_cut(lm[i])[0].to_highest_weight()[0] == rc_to_match:
                                res2.add(((*vlist,pw), vpl2_rc))
                
                
                    #print("After")
                    # vpl_new = vpl_new.vertical_cut(1)[1]
                    # if len(vpl_new) >= len(vpl):
                    #     vpl_new = vpl_new.vertical_cut(len(vpl) - 1)[0]
                    #vpl_new = vpl_new.rowrange(1)
                    #pretty_print(vpl_new)
                
                
                
                
                    # print("After")
                    # pretty_print(vpl_new)
                # vpl_new = vpl_new.extend(1)
                # print(f"After all reflections {shift_down=}")
                
                # last_code = 12000000
                # print("Pre cut")
                # pretty_print(vpl)
                # #pretty_print(vpl_new)
                
                # print("After cut")
                # pretty_print(vpl_new)
                #assert len(vpl_new) == len(vpl.perm.trimcode) and vpl_new.inv == vpl.inv

                
                
                # pw = tuple([a for a in vpl_new[0]])
                
                
                

                    
                    #     for i in range(len(vpl_new) - 1, 0, -1):
                    #         vpl_new_test = vpl_new.raising_operator(i)
                    #         if vpl_new_test is not None:                                
                    #             vpl_new = vpl_new_test
                    
                    
                    

                # v = vpl.perm
                #pw = tuple([v[ii] for ii in range(lm[i], len(v)) if ((ii > len(vpl_new.perm) and v[ii] == ii) or (ii <= len(vpl_new.perm) and v[ii] == vpl_new.perm[ii - 1]))])
                # if vpl_new.perm not in {vv[-1] for vv in vl}:
                #     print(f"OH NO NOT NULL OUT {lm[i] + 1}")
                #     pretty_print(vpl)
                #     print("ohas")
                #     pretty_print(vpl_new)
                #     print(f"{vl=}")
                #     raise AssertionError
                # pw = tuple(next(vv[0] for vv in vl if vv[-1] == vpl_new.perm))
                # #print(f"{vl=} {vpl_new.perm=} {pw=}")
                # # if len(pw) + v.inv != vpl.perm.inv:
                # #     continue
                # res2.add(((*vlist, pw), vpl_new))
                    

                # logger.debug(f"{vl=}")
                # rc_to_match = vpl.vertical_cut(lm[i])[0]
                # for pw, vpl2 in vl:
                #     for vpl2_rc in RCGraph.all_rc_graphs(vpl2, len(vpl)):
                #         if lm[i] + 1 >= len(vpl) or vpl2_rc.rowrange(lm[i] + 1) == vpl.rowrange(lm[i] + 1):
                #             if vpl2_rc.vertical_cut(lm[i])[0] == rc_to_match:
                #                 res2.add(((*vlist,pw), vpl2_rc))
                            
        res = res2
    if len(lm) == len(cn1w):
        return res
    res2 = set()
    for vlist, v_rc_0 in res:
        vp = v_rc_0
        vpl_list = divdiffable_rc(vp, c)
        if len(vpl_list) == 0:
            continue
        for vpl in vpl_list:
            res2.add((tuple(vlist), vpl))
    # logger.debug(f"{res2=}")
    return res2


def worker(nn, shared_recording_dict, lock, task_queue):
    import copy
    import sys
    from functools import cache
    from itertools import zip_longest

    import sympy
    from sympy import pretty_print

    from schubmult import CrystalGraphTensor, DSx, FreeAlgebra, Permutation, RCGraph, RCGraphRing, RootTableau, SchubertBasis, Sx


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
                            crystals[w_rc].add(tensor_lw)
                    
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
        from sympy import S
        try:
            (p, k, v, n) = task_queue.get(timeout=2)
        except Exception:
            break  # queue empty, exit
        with lock:
            if (p, k, v, n) in shared_recording_dict:
                pass
                if shared_recording_dict[(p, k, v, n)] is True:
                    #print(f"{(p, k, v, n)} already verified, returning.")
                    continue
                print(f"Previous failure on {(p, k, v, n)}, will retry.")
        
        
        rc_w_coprods = {}
        good = False
        #if True:
        #hw_tab = RCGraph.principal_rc(u, n - 1).to_highest_weight()[0]
        
        # mdom = Permutation.w0(n)#u.minimal_dominant_above()
        # w0 = mdom
        # # left diff
        # diff_perm = u * (~mdom)
        # poly_cache = {}
        # w0_prin = RCGraph.principal_rc(mdom, n)
        # sm = Sx.zero
        rc_ring = RCGraphRing()
        # W0 IS SPECIAL. THIS IS THE RC/CRYSTAL LEVEL DECOMPOSITION, NO ELEMENT
        # OTHER THAN w0 WORKS
        from itertools import zip_longest

        from schubmult import GeneratingSet, RCGraphRing, uncode
        rc_ring = RCGraphRing()

        good = True

        # lst = divdiffable_rc(w0_prin, ~v * mdom)
        # SKIP DDOM
        # if v == v.minimal_dominant_above():
        #     continue
        #bob = DSx(w0) * DSx(v, "z")
        # sm9 = Sx(w0) * Sx(v)
        # for boing, coeff in bob.items():
        #     # if coeff.expand() == S.Zero:
        #     #     continue
        #     print(f"w={boing} mu={mdom} coeff={coeff}")
        #     for v_rc in RCGraph.all_hw_rcs(v, n):
        #         print("rc=")
        #         pretty_print(v_rc)
        #         lst = dualpieri(mdom, v_rc=v_rc, w=boing)
        #     #print(coeff)
        #         print("result=")
        #         for bobb in lst:
        #             print(bobb[:-1])
        #             pretty_print(bobb[-1])
        #         #print([tuple(bobb[-1]) for bobb in lst])

        # if False:
        #     crystals = {}
        #     crystals[u] = {}
        #     for u_rc in RCGraph.all_rc_graphs(u, n):
        #         crystals[u][u_rc] = decompose_tensor_product(w0_prin, u_rc, n + 1)

        #     if u != v:
        #         crystals[v] = {}
        #         for v_rc in RCGraph.all_rc_graphs(v, n):
        #             crystals[v][v_rc] = decompose_tensor_product(w0_prin, v_rc, n + 1)
        # product = Sx(u) * Sx(v)
        # sm = Sx.zero
        # w0 = u.minimal_dominant_above()
        # pord = Sx(u) * Sx(v)
        # bongs = set()
        # bings = set()
        # sm0 = S.Zero
        # y = GeneratingSet("y")
        # z = GeneratingSet("z")
        # #check_elem = DSx([]).ring.from_dict({k: v for k, v in (DSx(u, "y") * DSx(v, "z")).items() if v.expand() != S.Zero})
        # # check_elem = DSx(u) * DSx(v, "z")
        # # x = Sx.genset
        # check_elem = DSx(u) * DSx(v, "z")
        # sm0 = DSx([]).ring.zero
        #for w in check_elem:
        # result = rc_ring.zero
        # complete_rc = RCGraph.principal_rc(uncode([*([0] * (k-1)), p]), n)
        def kk_insert_rc(rc, completerc):
            from schubmult import RCGraph
            rows = []
            if completerc.inv == 0:
                return rc
            for i, a in enumerate(completerc.length_vector):
                if a > 0:
                    rows.extend([i+1] * a)
            rows.reverse()
            return rc.kogan_kumar_insert(max(complete_rc.perm.descents()) + 1, rows)

        def kk_dist_mul(rc_elem, completerc):
            result = rc_elem.ring.zero
            for rc, coeff in rc_elem.items():
                new_rc = kk_insert_rc(rc, completerc)
                result += coeff * rc_elem.ring(new_rc)
            return result
        good = True
        for v_rc in RCGraph.all_lw_rcs(v, n):
            for complete_rc in RCGraph.all_rc_graphs(uncode([*([0] * (k-1)), p]), n):
                
                
                full_rc = kk_insert_rc(v_rc, complete_rc)
                for i in range(1, k):
                    elem1 = divdiff_rc_ring(rc_ring(full_rc), Permutation.ref_product(i))
                    #elem2_t = divdiff_rc_tensor_ring((rc_ring@rc_ring)(the_tensor.factors), Permutation.ref_product(i))
                    elem2 = rc_ring.zero
                    try:
                        v_rc_d, row = v_rc.exchange_property(i, return_row=True)
                        if row == i:
                            the_tensor = CrystalGraphTensor(v_rc_d, complete_rc.weight_reflection(i))
                            if the_tensor.raising_operator(i) is None:
                                while the_tensor is not None:
                                    elem2 += rc_ring(kk_insert_rc(*the_tensor.factors).normalize()  )
                                    the_tensor = the_tensor.lowering_operator(i)
                    except Exception:
                        pass
                    try:
                        assert all(val == 0 for val in (elem1 - elem2).values()), f"KK Dist fail {(p, k, v, n)} at reflection {i}"
                        print(f"{i=} {p=} {k=} {v_rc=} good")
                    except AssertionError as e:
                        print(e)
                        print("elem1")
                        pretty_print(elem1)
                        print("elem2")
                        pretty_print(elem2)
                        print("elem1 - elem2")
                        pretty_print(elem1 - elem2)
                        good = False

                    
            #dualpocket = dualpieri(u, v_rc,  w)
            # if len(dualpocket) > 0:
            #     #print(f"{u=} {w=} {v_rc=} {dualpocket=}")   
            #     for vlist, rc in dualpocket:
            #         toadd = S.One
            #         for i in range(len(vlist)):
            #             for j in range(len(vlist[i])):
            #                 toadd *= y[i + 1] - z[vlist[i][j]]
            #         sm0 += toadd * rc.polyvalue(y[len(vlist):], z) * DSx(w)
            #result += complete_mul(rc_ring(v_rc), p, k)
        # good = True
        
        # # diff = check_elem - sm0
        # # diff = DSx([]).ring.from_dict({k: sympy.sympify(vv).expand() for k, vv in diff.items() if sympy.sympify(vv).expand() != S.Zero})
        # #diff2 = check_elem - sm
        # try:
        #     #assert diff2 == Sx.zero, "Noit zor1"
        #     # from sympy import expand
        #     # assert all(expand(v99) == S.Zero for v99 in diff.values()), "Noit zor0"
        #     #assert diff == DSx([]).ring.zero

        #     #assert good
        #     #print(f"Coprod {rc.perm.trimcode}")
        #     # pretty_print(rc)
        #     # pretty_print(val)
        # #print("At least one success")
        #     #pretty_print(result)
        #     good = True
        # except AssertionError as e:
        #     # print(f"A fail {e=} {sm0=} {u=} {v=}")
        #     # print(f"{sm0=}")
        #     # print(f"{check_elem=}")
        #     # #print(f"{expand(sm0 - check_elem)=}")
        #     # # print(f"{u_rc.perm=}")
        #     # # print(f"{v=}")
        #     # print("diff=")
        #     # for w11, v0 in check_elem.items():
                
        #     #     if diff.get(w11, S.Zero).expand() == S.Zero:
        #     #         continue
        #     #     print("======")
        #     #     print("Actual")
        #     #     print(f"{w11}: {v0.expand()}")
        #     #     print("vs")
        #     #     print("Nope")
        #     #     print(f"{w11}: {sympy.sympify(sm0.get(w11, 0)).expand()}")
        #     # print(f"{sm0=}")
        #     # print(f"{check_elem=}")
        #     good = False
            
        #assert good, f"COMPLETE FAIL {w=}"
        if good:
            print(f"Success {(p, k, v, n)} at ", time.ctime())
            with lock:
                shared_recording_dict[(p, k, v, n)] = True
        else:
            with lock:
                shared_recording_dict[(p, k, v, n)] = False
            print(f"FAIL {(p, k, v, n)} at ", time.ctime())
    
        

def is_decomposable(w):
    for i in range(1, len(w) - 1):
        coset, w_J = w.coset_decomp(*list(range(1, i + 1)),*list(range(i + 2, len(w))))
        if coset.inv == 0 and set(w_J.code[:i+1]) != {0} and set(w_J.code[i+2:]) != {0}:
            return True
    return False

def main():
    from schubmult import DSx, Permutation, RCGraph, RootTableau, Sx, uncode

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
        dominant_only = True
        w0_only  = True
        sep_descs = False
        indec = False
        w0 = Permutation.w0(n)
        # for hw_tab in perms:
        #     if indec and is_decomposable(hw_tab):
        #         continue
        #     if (not dominant_only or hw_tab.minimal_dominant_above() == hw_tab) and (not w0_only or hw_tab == w0):
        #             # if indec and is_decomposable(perm):
        #             #     continue
        #             # if sep_descs:
        #             #     if hw_tab.inv == 0 or perm.inv == 0 or max(hw_tab.descents()) <= min(perm.descents()):
        #             #         task_queue.put((hw_tab, perm, n))
        #             # else:
        for k in range(n):
            for p in range(1, 2):
                for perm in perms:
                    task_queue.put((p, k, perm, n))

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

