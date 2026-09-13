# distutils: language = c++
# cython: language_level=3
"""schubmult_cpp: Python bindings for the C++ Schubert multiplication kernels in cpp/.

Coefficients are exchanged as symengine objects (the same RCP<const Basic> the Python
symengine wrapper holds), so nothing is expanded or re-parsed at the boundary.
"""

from libcpp.pair cimport pair
from libcpp.string cimport string
from libcpp.vector cimport vector

from symengine cimport rcp_const_basic
from symengine.lib.symengine_wrapper cimport Basic, c2py

import symengine

from schubmult.combinatorics.permutation import Permutation

cdef extern from "schubmult_api.h":
    ctypedef vector[int] PyPerm
    ctypedef vector[pair[PyPerm, int]] PyIntDict
    ctypedef vector[pair[PyPerm, rcp_const_basic]] PyExprDict
    PyIntDict api_schubmult_py(const PyIntDict& d, const PyPerm& v) except +
    ctypedef vector[rcp_const_basic] PyVars
    PyExprDict api_schubmult_double(const PyExprDict& d, const PyPerm& v, const PyVars& y, const PyVars& z) except +
    PyExprDict api_schubmult_q(const PyIntDict& d, const PyPerm& v, const PyVars& q) except +
    PyExprDict api_schubmult_q_double(const PyExprDict& d, const PyPerm& v, const PyVars& y, const PyVars& z, const PyVars& q) except +
    ctypedef rcp_const_basic (*ElemFuncCB)(void* ctx, int p, int k, const PyVars& y, const PyVars& z)
    PyExprDict api_schubmult_double_from_elems(const PyExprDict& d, const PyPerm& v, const PyVars& y, const PyVars& z, ElemFuncCB cb, void* ctx) except +
    PyExprDict api_schubmult_double_alt_from_elems(const PyExprDict& d, const PyPerm& v, const PyVars& y, const PyVars& z, ElemFuncCB cb, void* ctx) except +

cdef extern from "schub_common.h":
    int MAXN


# elem_func callback plumbing: the C++ kernel asks Python to build each distinct e_p(y | z)
# (memoized on the C++ side); an exception is parked here and re-raised after the kernel unwinds.
_pending_exc = []
_elem_funcs = {}  # id -> callable, kept alive while a kernel that may use it exists

cdef rcp_const_basic _elem_func_cb(void* ctx, int p, int k, const PyVars& y, const PyVars& z) noexcept with gil:
    cdef rcp_const_basic null
    cdef Basic b
    try:
        f = _elem_funcs[<size_t>ctx]
        yl = [c2py(y[i]) for i in range(y.size())]
        zl = [c2py(z[i]) for i in range(z.size())]
        b = symengine.sympify(f(p, k, yl, zl))
        return b.thisptr
    except BaseException as e:
        _pending_exc.append(e)
        return null


cdef object _with_elem_func(elem_func, run):
    """Register elem_func for the duration of run(ctx) and surface any exception it raised."""
    key = id(elem_func)
    _elem_funcs[key] = elem_func  # never evicted: the memo cache holding results keyed by it persists
    del _pending_exc[:]
    try:
        return run(key)
    except RuntimeError:
        if _pending_exc:
            raise _pending_exc[0]
        raise


cdef PyPerm _perm_in(perm) except *:
    cdef PyPerm p
    for x in perm:
        p.push_back(int(x))
    return p


cdef object _perm_out(const PyPerm& p):
    return Permutation(tuple(p))


cdef PyIntDict _int_dict_in(perm_dict) except *:
    cdef PyIntDict d
    cdef pair[PyPerm, int] e
    for k, v in perm_dict.items():
        e.first = _perm_in(k)
        e.second = int(v)  # raises TypeError for symbolic coefficients
        d.push_back(e)
    return d


cdef PyExprDict _expr_dict_in(perm_dict) except *:
    cdef PyExprDict d
    cdef pair[PyPerm, rcp_const_basic] e
    cdef Basic b
    for k, v in perm_dict.items():
        b = symengine.sympify(v)
        e.first = _perm_in(k)
        e.second = b.thisptr
        d.push_back(e)
    return d


cdef object _expr_dict_out(const PyExprDict& d):
    out = {}
    for i in range(d.size()):
        out[_perm_out(d[i].first)] = c2py(d[i].second)
    return out


cdef PyVars _vars_in(genset) except *:
    """The first MAXN+1 variables of a generating set (index i = i-th variable) as symengine objects.
    Entries past the end of the set are left null; the kernel errors only if it needs one."""
    cdef PyVars out
    cdef Basic b
    cdef rcp_const_basic null
    # ZeroGeneratingSet and friends have no length but index to 0 everywhere
    try:
        n = len(genset)
        if not isinstance(n, int):
            n = MAXN + 1
    except TypeError:
        n = MAXN + 1
    n = min(n, MAXN + 1)
    for i in range(n):
        b = symengine.sympify(genset[i])
        out.push_back(b.thisptr)
    for i in range(n, MAXN + 1):
        out.push_back(null)
    return out


def schubmult_py(perm_dict, v):
    """{Permutation: int} -> {Permutation: int}, the ordinary Schubert product."""
    cdef PyIntDict out = api_schubmult_py(_int_dict_in(perm_dict), _perm_in(v))
    return {_perm_out(out[i].first): out[i].second for i in range(out.size())}


def schubmult_double(perm_dict, v, var2, var3):
    """{Permutation: expr} -> {Permutation: expr}; var2 / var3 are generating sets (indexable, with len)."""
    return _expr_dict_out(api_schubmult_double(_expr_dict_in(perm_dict), _perm_in(v), _vars_in(var2), _vars_in(var3)))


def schubmult_q_fast(perm_dict, v, q_var):
    """Quantum product; input coefficients must be integers (the kernel is linear, so scale afterwards)."""
    return _expr_dict_out(api_schubmult_q(_int_dict_in(perm_dict), _perm_in(v), _vars_in(q_var)))


def schubmult_q_double_fast(perm_dict, v, var2, var3, q_var):
    return _expr_dict_out(api_schubmult_q_double(_expr_dict_in(perm_dict), _perm_in(v), _vars_in(var2), _vars_in(var3), _vars_in(q_var)))


def schubmult_double_from_elems(perm_dict, v, var2, var3, elem_func):
    """schubmult_double with elem_func(p, k, yvars, zvars) building the elementary symmetric objects
    (e.g. FactorialElemSym)."""
    cdef PyExprDict d = _expr_dict_in(perm_dict)
    cdef PyPerm vp = _perm_in(v)
    cdef PyVars y = _vars_in(var2)
    cdef PyVars z = _vars_in(var3)

    def run(size_t key):
        return _expr_dict_out(api_schubmult_double_from_elems(d, vp, y, z, _elem_func_cb, <void*>key))

    return _with_elem_func(elem_func, run)


def schubmult_double_alt_from_elems(perm_dict, v, var2, var3, elem_func):
    """schubmult_double_alt_from_elems_backwards: the pull-out-variable recursion with elem_func."""
    cdef PyExprDict d = _expr_dict_in(perm_dict)
    cdef PyPerm vp = _perm_in(v)
    cdef PyVars y = _vars_in(var2)
    cdef PyVars z = _vars_in(var3)

    def run(size_t key):
        return _expr_dict_out(api_schubmult_double_alt_from_elems(d, vp, y, z, _elem_func_cb, <void*>key))

    return _with_elem_func(elem_func, run)
