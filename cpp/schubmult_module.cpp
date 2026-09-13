// schubmult_module.cpp: the `schubmult.schubmult_cpp` extension module, written against the
// plain Python C API. Coefficients are Python symengine objects throughout (expr.h with
// SCHUB_PYEXPR), so building this needs only a C++17 compiler and Python headers.
//
//   schubmult_py(perm_dict, v)                          {Permutation: int} -> {Permutation: int}
//   schubmult_double(perm_dict, v, var2, var3[, elem_func])
//   schubmult_double_alt_from_elems(perm_dict, v, var2, var3, elem_func)
//   schubmult_q_fast(perm_dict, v, q_var)               integer input coefficients
//   schubmult_q_double_fast(perm_dict, v, var2, var3, q_var)
//
// Permutations are anything iterable over ints (one-line notation); results are keyed by
// schubmult.combinatorics.permutation.Permutation. A RuntimeError signals "exceeds MAXN".

#define SCHUB_PYEXPR
#include "schubmult_api.h"

namespace {

// Lazily imported: the module lives inside the schubmult package, so importing Permutation at
// module-init time would be circular.
PyObject* permutation_class() {
    static PyObject* cls = nullptr;
    if (!cls) {
        Expr mod = pyexpr_detail::checked(PyImport_ImportModule("schubmult.combinatorics.permutation"));
        cls = pyexpr_detail::checked(PyObject_GetAttrString(mod.get(), "Permutation")).release();
    }
    return cls;
}

PyPerm perm_in(PyObject* obj) {
    Expr seq = pyexpr_detail::checked(PySequence_Fast(obj, "permutation must be iterable"));
    Py_ssize_t n = PySequence_Fast_GET_SIZE(seq.get());
    PyPerm p;
    p.reserve((size_t)n);
    for (Py_ssize_t i = 0; i < n; ++i) {
        long x = PyLong_AsLong(PySequence_Fast_GET_ITEM(seq.get(), i));
        if (x == -1 && PyErr_Occurred()) throw PyErrorSet();
        p.push_back((int)x);
    }
    return p;
}

Expr perm_out(const PyPerm& p) {
    Expr tup = pyexpr_detail::checked(PyTuple_New((Py_ssize_t)p.size()));
    for (size_t i = 0; i < p.size(); ++i) PyTuple_SET_ITEM(tup.get(), (Py_ssize_t)i, pyexpr_detail::checked(PyLong_FromLong(p[i])).release());
    return pyexpr_detail::checked(PyObject_CallFunctionObjArgs(permutation_class(), tup.get(), nullptr));
}

PyIntDict int_dict_in(PyObject* d) {
    PyIntDict out;
    PyObject *k, *v;
    Py_ssize_t pos = 0;
    while (PyDict_Next(d, &pos, &k, &v)) {
        Expr idx = pyexpr_detail::checked(PyNumber_Index(v));  // symbolic coefficients raise TypeError
        long c = PyLong_AsLong(idx.get());
        if (c == -1 && PyErr_Occurred()) throw PyErrorSet();
        out.push_back({perm_in(k), (int)c});
    }
    return out;
}

PyExprDict expr_dict_in(PyObject* d) {
    PyExprDict out;
    PyObject *k, *v;
    Py_ssize_t pos = 0;
    while (PyDict_Next(d, &pos, &k, &v)) out.push_back({perm_in(k), ex_sympify(Expr(v, false))});
    return out;
}

Expr int_dict_out(const PyIntDict& d) {
    Expr r = pyexpr_detail::checked(PyDict_New());
    for (const auto& kv : d) {
        Expr key = perm_out(kv.first), val = pyexpr_detail::checked(PyLong_FromLong(kv.second));
        if (PyDict_SetItem(r.get(), key.get(), val.get()) < 0) throw PyErrorSet();
    }
    return r;
}

Expr expr_dict_out(const PyExprDict& d) {
    Expr r = pyexpr_detail::checked(PyDict_New());
    for (const auto& kv : d) {
        Expr key = perm_out(kv.first);
        if (PyDict_SetItem(r.get(), key.get(), kv.second.get()) < 0) throw PyErrorSet();
    }
    return r;
}

// The first MAXN+1 variables of a generating set (index i = i-th variable); entries past its end
// are left null and only error if the kernel needs them. Sets without a usable len() (e.g.
// ZeroGeneratingSet) are indexed all the way.
PyVars vars_in(PyObject* genset) {
    PyVars out;
    Py_ssize_t n = MAXN + 1;
    Py_ssize_t len = PyObject_Length(genset);
    if (len < 0)
        PyErr_Clear();
    else
        n = std::min(n, len);
    for (Py_ssize_t i = 0; i < n; ++i) {
        Expr item = pyexpr_detail::checked(PySequence_GetItem(genset, i));
        out.push_back(ex_sympify(item));
    }
    out.resize(MAXN + 1);
    return out;
}

// Runs f, converting C++ failures to Python exceptions.
template <typename F>
PyObject* guarded(F f) {
    try {
        return f().release();
    } catch (const PyErrorSet&) {
        return nullptr;  // the Python error is already set
    } catch (const std::exception& e) {
        if (!PyErr_Occurred()) PyErr_SetString(PyExc_RuntimeError, e.what());
        return nullptr;
    }
}

PyObject* py_schubmult_py(PyObject*, PyObject* args) {
    PyObject *d, *v;
    if (!PyArg_ParseTuple(args, "O!O", &PyDict_Type, &d, &v)) return nullptr;
    return guarded([&] { return int_dict_out(api_schubmult_py(int_dict_in(d), perm_in(v))); });
}

PyObject* py_schubmult_double(PyObject*, PyObject* args) {
    PyObject *d, *v, *var2, *var3, *elem_func = Py_None;
    if (!PyArg_ParseTuple(args, "O!OOO|O", &PyDict_Type, &d, &v, &var2, &var3, &elem_func)) return nullptr;
    return guarded([&] {
        Expr ef = elem_func == Py_None ? Expr() : Expr(elem_func, false);
        return expr_dict_out(api_schubmult_double(expr_dict_in(d), perm_in(v), vars_in(var2), vars_in(var3), ef));
    });
}

PyObject* py_schubmult_double_alt_from_elems(PyObject*, PyObject* args) {
    PyObject *d, *v, *var2, *var3, *elem_func;
    if (!PyArg_ParseTuple(args, "O!OOOO", &PyDict_Type, &d, &v, &var2, &var3, &elem_func)) return nullptr;
    return guarded([&] { return expr_dict_out(api_schubmult_double_alt_from_elems(expr_dict_in(d), perm_in(v), vars_in(var2), vars_in(var3), Expr(elem_func, false))); });
}

PyObject* py_schubmult_q_fast(PyObject*, PyObject* args) {
    PyObject *d, *v, *q;
    if (!PyArg_ParseTuple(args, "O!OO", &PyDict_Type, &d, &v, &q)) return nullptr;
    return guarded([&] { return expr_dict_out(api_schubmult_q(int_dict_in(d), perm_in(v), vars_in(q))); });
}

PyObject* py_schubmult_q_double_fast(PyObject*, PyObject* args) {
    PyObject *d, *v, *var2, *var3, *q;
    if (!PyArg_ParseTuple(args, "O!OOOO", &PyDict_Type, &d, &v, &var2, &var3, &q)) return nullptr;
    return guarded([&] { return expr_dict_out(api_schubmult_q_double(expr_dict_in(d), perm_in(v), vars_in(var2), vars_in(var3), vars_in(q))); });
}

PyMethodDef methods[] = {
    {"schubmult_py", py_schubmult_py, METH_VARARGS, "Ordinary Schubert product: {Permutation: int} -> {Permutation: int}."},
    {"schubmult_double", py_schubmult_double, METH_VARARGS, "Double Schubert product; var2/var3 are generating sets, optional elem_func builds the elementary symmetric objects."},
    {"schubmult_double_alt_from_elems", py_schubmult_double_alt_from_elems, METH_VARARGS, "Pull-out-variable recursion with elem_func (positivity)."},
    {"schubmult_q_fast", py_schubmult_q_fast, METH_VARARGS, "Quantum Schubert product; integer input coefficients."},
    {"schubmult_q_double_fast", py_schubmult_q_double_fast, METH_VARARGS, "Quantum double Schubert product."},
    {nullptr, nullptr, 0, nullptr},
};

PyModuleDef moduledef = {PyModuleDef_HEAD_INIT, "schubmult_cpp", "C++ Schubert multiplication kernels.", -1, methods, nullptr, nullptr, nullptr, nullptr};

}  // namespace

PyMODINIT_FUNC PyInit_schubmult_cpp() {
    PyObject* m = PyModule_Create(&moduledef);
    if (m) PyModule_AddIntConstant(m, "MAXN", MAXN);
    return m;
}
