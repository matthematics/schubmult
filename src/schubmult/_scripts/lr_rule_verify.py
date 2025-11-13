import argparse
import base64
import pickle
import shutil
import sys
import time
import os
import json
from multiprocessing import Event, Manager, Process


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
        rc_ring = RCGraphRing()
        tring = rc_ring @ rc_ring

        crystals = {}

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

        # cut_left_rc = left_rc.vertical_cut(n - 2)[0]
        # cut_u = u_rc.vertical_cut(n - 2)[0]

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
        # used_tensors = set()
        for rc_w_cut, tensor_elems in cut_crystals.items():
            up_rc = rc_ring(RCGraph.one_row(top_row_left + top_row_right)) * rc_ring(rc_w_cut)
            for w_rc, coeff in up_rc.items():
                assert coeff == 1
                high_weight = w_rc.to_highest_weight()[0].crystal_weight
                low_weight = w_rc.to_lowest_weight()[0].crystal_weight
                for t_elem in tensor_elems:
                    # up_rc =  rc_ring(rc_w_cut) * rc_ring(RCGraph.one_row(len(left_rc[-1]) + len(u_rc[-1])))
                    up_tensor = tring((RCGraph.one_row(top_row_left), RCGraph.one_row(top_row_right))) * tring(t_elem.factors)
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
                        # tensor_dom = CrystalGraphTensor(dom_rc, u_rc2)
                        tensor_hw = tensor.to_highest_weight()[0]
                        tensor_lw = tensor.to_lowest_weight()[0]
                        # if tensor_hw in used_tensors:
                        #     continue
                        # low_tensor_weight = tuple([a + b for a,b in zip(left_rc.to_lowest_weight()[0].length_vector, tensor_lw.factors[1].length_vector)])
                        low_tensor_weight = tensor_lw.crystal_weight

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
                            # used_tensors.add(tensor_hw)
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

    while True:
        from sympy import S

        try:
            item = task_queue.get()
            if item is None:
                break
            (key, check_val) = item
            (u, v, w, n) = key
            check_val = Sx.from_dict(check_val)
        except Exception:
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

        mdom = Permutation.w0(n)  # u.minimal_dominant_above()

        w0_prin = RCGraph.principal_rc(mdom, n)

        # rc_ring = RCGraphRing()

        good = True

        sm0 = S.Zero
    
        

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
        good = True
        #check_elem = Sx(u) * Sx(v)
        for w in check_val:
            if not good:
                break
            for tensor_hw in hw_tensors:
                if not good:
                    break
                for inner_tensor in tensor_hw.full_crystal:
                    if not good:
                        break
                    v_rc = inner_tensor.factors[1]
                    dualpocket = v_rc.dualpieri(u, w)
                    if len(dualpocket) > 0:
                        if tensor_hw in visited:
                            good = False
                            break
                        # if len() == 0:
                        #     good = False
                        #     break
                        rc_high = RCGraph.all_rc_graphs(w, n, weight=tensor_hw.crystal_weight)
                        rc_low = RCGraph.all_rc_graphs(w, n, weight=inner_tensor.to_lowest_weight()[0].crystal_weight)
                        good = False
                        for rc_lw in rc_low:
                            if not rc_lw.is_lowest_weight or rc_lw.to_highest_weight()[0] not in rc_high:
                                continue
                            good = True
                        if not good:
                            break
                        visited.add(tensor_hw)
                        dualps[w] = dualps.get(w, set())
                        for vlist, perm_list, rc in dualpocket:
                            # if (vlist, perm_list, rc) in dualps[w]:
                            #     continue

                            dualps[w].add((vlist, perm_list, rc))
                            toadd = S.One
                            for i in range(len(vlist)):
                                for j in range(len(vlist[i])):
                                    toadd *= y[i + 1] - z[vlist[i][j]]
                            sm0 += Sx(w)#vlist.polyvalue(y, z) * rc.polyvalue(y[len(vlist) :], z)# * DSx(w)

        ########################333
        #check_elem = check_elem.ring.from_dict({u: check_elem[u]})
        #############3
        from sympy import expand
        diff = expand(check_val - sm0)
        # diff = DSx([]).ring.from_dict({k: sympy.sympify(vv).expand() for k, vv in diff.items() if sympy.sympify(vv).expand() != S.Zero})

        try:
            

            # assert all(expand(v99) == S.Zero for v99 in diff.values()), "Noit zor0"
            assert good
            assert diff == S.Zero, "Noit zor0"
            #assert len(visited) == len(hw_tensors), "Missed tensors"

            
            
        except AssertionError as e:
            print(f"A fail {e=} {sm0=} {u=} {v=} {w=}")
            print(f"{sm0=}")
            print(f"{diff=}")
            # for w11, v0 in check_elem.items():
            #     if diff.get(w11, S.Zero).expand() == S.Zero:
            #         continue
                # print("======")
                # print("Actual")
                # print(f"{w11}: {v0.expand()}")
                # print("vs")
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
        "irreducible": True,
        "skip_id": True,
        "molevsagan": False,
    }

    bool_group("dominant_only", None)
    bool_group("w0_only", None)
    bool_group("base_level", None)
    bool_group("sep_descs", None)
    bool_group("irreducible", None)
    bool_group("skip_id", None)
    bool_group("molevsagan", None)

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
        if dominant_only:
            print("Only verifying dominant multipliers.")
        if w0_only:
            print("Only verifying w0 RC.")
        if sep_descs:
            print("Testing only separated descents.")
        if skip_id:
            print("Skipping identity.")

        # queue population (unchanged logic), but skip enqueuing items already True in shared_recording_dict
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
                            key = (hw_tab, perm, n)
                            if shared_recording_dict.get(key) is True:
                                continue
                            task_queue.put((hw_tab, perm, n))
                    else:
                        if molevsagan:
                            check_elem = DSx(hw_tab) * DSx(perm, "z")
                        else:
                            check_elem = Sx(hw_tab) * Sx(perm)
                        #keys = set(check_elem.keys())
                        #for w in keys:
                        val = check_elem
                        #key = (hw_tab, perm, w, n)
                        key = (hw_tab, perm, None, n)
                        if shared_recording_dict.get(key) is True:
                            continue
                        task_queue.put((key, dict(val)))
                        # for w in keys:
                        #     del check_elem[w]

        # Start fixed number of workers and place sentinels so they exit cleanly
        from schubmult import y, z
        workers = []
        for _ in range(num_processors):
            task_queue.put(None)  # sentinel
            
            p = Process(target=worker, args=(n, shared_recording_dict, lock, task_queue, y, z))
            p.start()
            workers.append(p)
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
