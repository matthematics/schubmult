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

    # try:
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    for perm in perms:
        mod = try_lr_module(perm)
        word_of_perm = []
        weight_of_perm = []
        for i, v in enumerate(perm.trimcode):
            word_of_perm.extend(list(range(i + v, i, -1)))
            weight_of_perm.extend([i + 1] * v)
        for (rc1, rc2), coeff in mod.items():
            print(f"{rc1.perm, rc2.perm} for {perm} with coeff {coeff}")
            print(f"{[*rc1.perm_word(),*rc2.perm_word()]}")
            print(f"{word_of_perm}")
            print("and")
            print(f"{[*rc1.weight_word(),*rc2.weight_word()]}")
            print(f"{weight_of_perm}")
            print("----")

if __name__ == "__main__":
    main()
