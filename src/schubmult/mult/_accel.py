"""C++ multiplication kernels (the ``schubmult_cpp`` extension, built from ``cpp/`` by setup.py).

``schubmult_py``, ``schubmult_double``, ``schubmult_q_fast``, ``schubmult_q_double_fast`` and the
``*_from_elems`` kernels dispatch here. The extension is a required part of the package; the
pure-Python kernels remain only as the fallback for permutations beyond the compiled MAXN (the
wrappers return ``None`` in that case). Set ``SCHUBMULT_NO_CPP=1`` to force the Python kernels.
"""

import os

if os.environ.get("SCHUBMULT_NO_CPP"):
    _cpp = None
else:
    try:
        from schubmult import schubmult_cpp as _cpp
    except ImportError as e:
        raise ImportError("the schubmult_cpp extension is missing; reinstall the package (`pip install .`), which compiles it") from e

available = _cpp is not None


def _is_int(x):
    """Whether ``x`` is (numerically) an integer."""
    try:
        return int(x) == x
    except Exception:
        return False


def _call(f, *args):
    """Call the extension function ``f``; return ``None`` on its ``RuntimeError`` (exceeds ``MAXN``) signal."""
    # RuntimeError is the extension's "exceeds MAXN" signal; let the caller fall back
    try:
        return f(*args)
    except RuntimeError:
        return None


def schubmult_py(perm_dict, v):
    """C++ ``schubmult_py``; returns ``None`` if any coefficient is non-integer or the size exceeds ``MAXN``."""
    if not all(_is_int(c) for c in perm_dict.values()):
        return None
    return _call(_cpp.schubmult_py, {k: int(c) for k, c in perm_dict.items()}, v)


def schubmult_double(perm_dict, v, var2, var3):
    """C++ ``schubmult_double``; returns ``None`` if the size exceeds ``MAXN``."""
    return _call(_cpp.schubmult_double, perm_dict, v, var2, var3)


def schubmult_q_fast(perm_dict, v, q_var):
    """C++ ``schubmult_q_fast``; non-integer coefficients are handled by linearity, one key at a time."""
    if all(_is_int(c) for c in perm_dict.values()):
        return _call(_cpp.schubmult_q_fast, {k: int(c) for k, c in perm_dict.items()}, v, q_var)
    # the kernel is linear in the coefficients
    ret = {}
    for k, c in perm_dict.items():
        part = _call(_cpp.schubmult_q_fast, {k: 1}, v, q_var)
        if part is None:
            return None
        for w, val in part.items():
            ret[w] = ret.get(w, 0) + c * val
    return ret


def schubmult_q_double_fast(perm_dict, v, var2, var3, q_var):
    """C++ ``schubmult_q_double_fast``; returns ``None`` if the size exceeds ``MAXN``."""
    return _call(_cpp.schubmult_q_double_fast, perm_dict, v, var2, var3, q_var)


def schubmult_double_from_elems(perm_dict, v, var2, var3, elem_func):
    """C++ ``schubmult_double`` with a custom elementary-symmetric ``elem_func``; ``None`` if size exceeds ``MAXN``."""
    return _call(_cpp.schubmult_double, perm_dict, v, var2, var3, elem_func)


def schubmult_double_alt_from_elems(perm_dict, v, var2, var3, elem_func):
    """C++ ``schubmult_double_alt_from_elems``; ``None`` if size exceeds ``MAXN``."""
    return _call(_cpp.schubmult_double_alt_from_elems, perm_dict, v, var2, var3, elem_func)
