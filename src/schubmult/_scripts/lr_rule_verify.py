# LR rule verification script

import base64
import json
import os
import pickle
import shutil
import sys
import time
from multiprocessing import Event, Manager, Process


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
    except Exception:
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
        except Exception:
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


def worker(nn, shared_recording_dict, lock, task_queue):  # noqa: ARG001
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
            (u, v, n) = task_queue.get(timeout=2)
        except Exception:
            break  # queue empty, exit
        with lock:
            if (u, v, n) in shared_recording_dict:
                if shared_recording_dict[(u, v, n)] is True:
                    # print(f"{(u, v, n)} already verified, returning.")
                    continue
                print(f"Previous failure on {(u, v, n)}, will retry.")

        good = False
        # if True:
        # hw_tab = RCGraph.principal_rc(u, n - 1).to_highest_weight()[0]

        mdom = Permutation.w0(n)  # u.minimal_dominant_above()

        w0_prin = RCGraph.principal_rc(mdom, n)

        from schubmult import GeneratingSet

        # rc_ring = RCGraphRing()

        good = True

        sm0 = S.Zero
        y = GeneratingSet("y")
        z = GeneratingSet("z")

        check_elem = DSx(u) * DSx(v, "z")

        sm0 = DSx([]).ring.zero
        dualps = {}
        for w in check_elem:
            for v_rc in RCGraph.all_rc_graphs(v, n):
                dualpocket = v_rc.dualpieri(u, w)
                if len(dualpocket) > 0:
                    dualps[w] = dualps.get(w, set())
                    for vlist, perm_list, rc in dualpocket:
                        if (vlist, perm_list, rc) in dualps[w]:
                            continue

                        dualps[w].add((vlist, perm_list, rc))
                        toadd = S.One
                        for i in range(len(vlist)):
                            for j in range(len(vlist[i])):
                                toadd *= y[i + 1] - z[vlist[i][j]]
                        sm0 += toadd * rc.polyvalue(y[len(vlist) :], z) * DSx(w)

        good = True

        diff = check_elem - sm0
        diff = DSx([]).ring.from_dict({k: sympy.sympify(vv).expand() for k, vv in diff.items() if sympy.sympify(vv).expand() != S.Zero})

        try:
            from sympy import expand

            assert all(expand(v99) == S.Zero for v99 in diff.values()), "Noit zor0"

            good = True
        except AssertionError as e:
            print(f"A fail {e=} {sm0=} {u=} {v=}")
            print(f"{sm0=}")
            print(f"{check_elem=}")
            print("diff=")
            for w11, v0 in check_elem.items():
                if diff.get(w11, S.Zero).expand() == S.Zero:
                    continue
                print("======")
                print("Actual")
                print(f"{w11}: {v0.expand()}")
                print("vs")
                for dual in dualps.get(w11, set()):
                    for domp in dual:
                        print(domp)
                print(f"{w11}: {sm0.get(w11)}")
            good = False

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
        coset, w_J = w.coset_decomp(*list(range(1, i + 1)), *list(range(i + 2, len(w))))
        if coset.inv == 0 and set(w_J.code[: i + 1]) != {0} and set(w_J.code[i + 2 :]) != {0}:
            return True
    return False


def main():
    from schubmult import DSx, Permutation, RCGraph, RootTableau, Sx, uncode  # noqa: F401

    try:
        n = int(sys.argv[1])
        filename = None
        verification_filename = None
        num_processors = int(sys.argv[2])
        if len(sys.argv) > 3:
            filename = sys.argv[2]
            verification_filename = filename + ".verification"
        else:
            print("No verification filename provided, not saving results to disk.", file=sys.stderr)

    except (IndexError, ValueError):
        print("Usage: lr_rule_verify n num_processors <filename>", file=sys.stderr)
        print("filename.verification.json is used for loading/saving verification results if provided", file=sys.stderr)
        sys.exit(1)

    perms = Permutation.all_permutations(n)

    perms.sort(key=lambda p: (p.inv, p.trimcode))

    with Manager() as manager:
        shared_recording_dict = manager.dict()
        lock = manager.Lock()
        stop_event = Event()

        # Load recording dict from JSON only
        if verification_filename is not None:
            loaded_recording = safe_load_recording(verification_filename)
            if loaded_recording:
                shared_recording_dict.update(loaded_recording)

            recording_saver_proc = Process(target=recording_saver, args=(shared_recording_dict, lock, verification_filename, stop_event))
            recording_saver_proc.start()

        task_queue = manager.Queue()
        dominant_only = True
        w0_only = False
        sep_descs = False
        indec = False
        w0 = Permutation.w0(n)
        for hw_tab in perms:
            if indec and is_decomposable(hw_tab):
                continue
            if (not dominant_only or hw_tab.minimal_dominant_above() == hw_tab) and (not w0_only or hw_tab == w0):
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
