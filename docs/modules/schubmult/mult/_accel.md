<a id="schubmult.mult._accel"></a>

# schubmult.mult.\_accel

C++ multiplication kernels (the ``schubmult_cpp`` extension, built from ``cpp/`` by setup.py).

``schubmult_py``, ``schubmult_double``, ``schubmult_q_fast``, ``schubmult_q_double_fast`` and the
``*_from_elems`` kernels dispatch here. The extension is a required part of the package; the
pure-Python kernels remain only as the fallback for permutations beyond the compiled MAXN (the
wrappers return ``None`` in that case). Set ``SCHUBMULT_NO_CPP=1`` to force the Python kernels.

<a id="schubmult.mult._accel.schubmult_py"></a>

#### schubmult\_py

```python
def schubmult_py(perm_dict, v)
```

C++ ``schubmult_py``; returns ``None`` if any coefficient is non-integer or the size exceeds ``MAXN``.

<a id="schubmult.mult._accel.schubmult_double"></a>

#### schubmult\_double

```python
def schubmult_double(perm_dict, v, var2, var3)
```

C++ ``schubmult_double``; returns ``None`` if the size exceeds ``MAXN``.

<a id="schubmult.mult._accel.schubmult_q_fast"></a>

#### schubmult\_q\_fast

```python
def schubmult_q_fast(perm_dict, v, q_var)
```

C++ ``schubmult_q_fast``; non-integer coefficients are handled by linearity, one key at a time.

<a id="schubmult.mult._accel.schubmult_q_double_fast"></a>

#### schubmult\_q\_double\_fast

```python
def schubmult_q_double_fast(perm_dict, v, var2, var3, q_var)
```

C++ ``schubmult_q_double_fast``; returns ``None`` if the size exceeds ``MAXN``.

<a id="schubmult.mult._accel.schubmult_double_from_elems"></a>

#### schubmult\_double\_from\_elems

```python
def schubmult_double_from_elems(perm_dict, v, var2, var3, elem_func)
```

C++ ``schubmult_double`` with a custom elementary-symmetric ``elem_func``; ``None`` if size exceeds ``MAXN``.

<a id="schubmult.mult._accel.schubmult_double_alt_from_elems"></a>

#### schubmult\_double\_alt\_from\_elems

```python
def schubmult_double_alt_from_elems(perm_dict, v, var2, var3, elem_func)
```

C++ ``schubmult_double_alt_from_elems``; ``None`` if size exceeds ``MAXN``.

