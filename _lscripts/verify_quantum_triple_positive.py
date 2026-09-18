"""Parallel, restartable verification of positive representations for quantum
double Schubert triple products.

Each individual case is keyed by ``(perm1, perm2, perm3)`` and logged as its
own JSON record (JSON Lines format) to a log file. There is no shared,
in-memory results dictionary: workers pull independent ``(perm1, perm2)``
tasks off a queue, compute the triple product, verify each ``perm3`` term,
and append one record per case directly to the log file. This keeps workers
fully independent and makes the log the single source of truth.

Restartability: on startup the log file (if it exists) is scanned once to
recover the set of ``(perm1, perm2, perm3)`` cases already verified, plus the
set of ``(perm1, perm2)`` pairs that were fully completed. Both are skipped
on subsequent runs.

Overwrite safety: the log file is only ever opened in append mode, so a
crashed or restarted run can never truncate or clobber previously recorded
results. The default log filename is derived from ``n`` so runs for
different ``n`` do not collide by default.
"""

import argparse
import itertools
import json
import multiprocessing as mp
import os
import sys
import time


def _perm_key(perm):
    """JSON-friendly representation of a Permutation (its array form)."""
    return list(perm)


def _load_progress(log_path):
    """
    Scan an existing JSONL log and recover:
      - completed_triples: set of (perm1, perm2, perm3) tuples already verified
      - completed_pairs: set of (perm1, perm2) tuples fully completed
    Corrupt/partial trailing lines (e.g. from a killed process) are ignored.
    """
    completed_triples = set()
    completed_pairs = set()
    if not os.path.exists(log_path):
        return completed_triples, completed_pairs

    with open(log_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            event = rec.get("event")
            if event == "case" and rec.get("verified") is True:
                key = (tuple(rec["perm1"]), tuple(rec["perm2"]), tuple(rec["perm3"]))
                completed_triples.add(key)
            elif event == "pair_complete":
                pkey = (tuple(rec["perm1"]), tuple(rec["perm2"]))
                completed_pairs.add(pkey)
    return completed_triples, completed_pairs


def _append_record(log_path, lock, record):
    line = json.dumps(record) + "\n"
    with lock:
        with open(log_path, "a") as f:
            f.write(line)
            f.flush()
            os.fsync(f.fileno())


def _worker(worker_id, task_queue, log_path, lock, completed_triples, n, stop_event, fail_event, stop_on_failure):
    from schubmult.abc import q, y, z
    from schubmult.mult.quantum_double import can_classically_reduce, schubmult_q_double

    var2, var3 = y, z

    while True:
        if stop_event.is_set():
            break
        task = task_queue.get()
        if task is None:
            break
        perm1, perm2 = task

        if stop_on_failure and fail_event.is_set():
            continue

        p1key = _perm_key(perm1)
        p2key = _perm_key(perm2)

        vdict = schubmult_q_double({perm1: 1}, perm2, var2, var3, q)
        pair_ok = True
        for perm3, val in vdict.items():
            p3key = _perm_key(perm3)
            triple_key = (tuple(p1key), tuple(p2key), tuple(p3key))
            if triple_key in completed_triples:
                continue

            verified = can_classically_reduce(perm1, perm2, perm3, val, q, filter_dominant=True) # bool(q_posify(perm1, perm2, perm3, val, var2, var3, q, False, verify_only=True))

            _append_record(
                log_path,
                lock,
                {
                    "event": "case",
                    "perm1": p1key,
                    "perm2": p2key,
                    "perm3": p3key,
                    "verified": verified,
                    "val": str(val),
                    "worker": worker_id,
                    "time": time.ctime(),
                },
            )

            if not verified:
                pair_ok = False
                fail_event.set()
                print(f"FAIL {perm1=} {perm2=} {perm3=} {val=}", flush=True)
                if stop_on_failure:
                    break

        if pair_ok:
            _append_record(
                log_path,
                lock,
                {
                    "event": "pair_complete",
                    "perm1": p1key,
                    "perm2": p2key,
                    "n": n,
                    "time": time.ctime(),
                },
            )
            print(f"Verified positive representation for {perm1=} {perm2=}", flush=True)


def _queue_producer(task_queue, perms, num_processors, completed_pairs):
    for perm1, perm2 in itertools.product(perms, repeat=2):
        pkey = (tuple(_perm_key(perm1)), tuple(_perm_key(perm2)))
        if pkey in completed_pairs:
            continue
        task_queue.put((perm1, perm2))
    for _ in range(num_processors):
        task_queue.put(None)


def main():
    parser = argparse.ArgumentParser(description="Parallel restartable verification of quantum triple positivity")
    parser.add_argument("n", type=int)
    parser.add_argument("num_processors", type=int)
    parser.add_argument("--log-file", dest="log_file", help="path to the JSONL log file (default derived from n)")
    parser.add_argument("--stop-on-failure", action="store_true", help="stop scheduling new work after the first failure")
    args = parser.parse_args()

    from schubmult import Permutation

    sys.setrecursionlimit(1000000)

    n = args.n
    num_processors = args.num_processors
    log_path = args.log_file or f"quantum_triple_positive_n{n}.verification.jsonl"

    completed_triples, completed_pairs = _load_progress(log_path)
    if completed_triples or completed_pairs:
        print(f"Resuming from {log_path}: {len(completed_pairs)} pairs and {len(completed_triples)} cases already verified.", flush=True)

    perms = Permutation.all_permutations(n)

    manager = mp.Manager()
    task_queue = manager.Queue()
    lock = manager.Lock()
    stop_event = manager.Event()
    fail_event = manager.Event()

    producer = mp.Process(target=_queue_producer, args=(task_queue, perms, num_processors, completed_pairs))
    producer.start()

    workers = []
    for worker_id in range(num_processors):
        p = mp.Process(
            target=_worker,
            args=(worker_id, task_queue, log_path, lock, completed_triples, n, stop_event, fail_event, args.stop_on_failure),
        )
        p.start()
        workers.append(p)

    producer.join()
    for p in workers:
        p.join()

    if fail_event.is_set():
        print("Verification FAILED for one or more cases; see", log_path, flush=True)
        sys.exit(1)

    print("All cases verified.", flush=True)


if __name__ == "__main__":
    main()

                