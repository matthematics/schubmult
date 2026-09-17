<a id="schubmult._scripts.verify_quantum_triple_positive"></a>

# schubmult.\_scripts.verify\_quantum\_triple\_positive

Parallel, restartable verification of positive representations for quantum
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

