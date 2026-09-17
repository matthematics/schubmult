<a id="schubmult.mult._accel"></a>

# schubmult.mult.\_accel

C++ multiplication kernels (the ``schubmult_cpp`` extension, built from ``cpp/`` by setup.py).

``schubmult_py``, ``schubmult_double``, ``schubmult_q_fast``, ``schubmult_q_double_fast`` and the
``*_from_elems`` kernels dispatch here. The extension is a required part of the package; the
pure-Python kernels remain only as the fallback for permutations beyond the compiled MAXN (the
wrappers return ``None`` in that case). Set ``SCHUBMULT_NO_CPP=1`` to force the Python kernels.

