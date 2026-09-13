// expr.h: the symbolic-coefficient type used by the double / quantum-double kernels, with two
// interchangeable backends selected at compile time:
//
//   default        SymEngine C++ (RCP<const Basic>) — for the standalone CLI tools.
//   SCHUB_PYEXPR   Python objects (PyObject*) driven through the Python C API — for the
//                  schubmult_cpp extension. Coefficients are then whatever the `symengine`
//                  Python package hands back, so the extension needs no SymEngine C++ headers
//                  or library and links against nothing but Python. Any Python exception is
//                  carried as a C++ exception (PyErrorSet) to the module boundary.
//
// The kernels only use the small ex_* vocabulary below.

#pragma once

#include <stdexcept>
#include <string>
#include <vector>

#ifdef SCHUB_PYEXPR
// ---------------------------------------------------------------------------
// Python-object backend
// ---------------------------------------------------------------------------
#define PY_SSIZE_T_CLEAN
#include <Python.h>

struct PyErrorSet : std::runtime_error {
    PyErrorSet() : std::runtime_error("Python exception") {}
};

// Owning reference to a Python object (null allowed).
class Expr {
    PyObject* p_ = nullptr;

public:
    Expr() = default;
    explicit Expr(PyObject* p, bool steal = true) : p_(p) {
        if (p_ && !steal) Py_INCREF(p_);
    }
    Expr(const Expr& o) : p_(o.p_) {
        if (p_) Py_INCREF(p_);
    }
    Expr(Expr&& o) noexcept : p_(o.p_) { o.p_ = nullptr; }
    Expr& operator=(Expr o) noexcept {
        std::swap(p_, o.p_);
        return *this;
    }
    ~Expr() {
        if (p_) Py_DECREF(p_);
    }
    bool is_null() const { return p_ == nullptr; }
    PyObject* get() const { return p_; }
    PyObject* release() {
        PyObject* r = p_;
        p_ = nullptr;
        return r;
    }
};

typedef std::vector<Expr> ExprVec;

namespace pyexpr_detail {
// Result of a Python call; throws if the call failed.
inline Expr checked(PyObject* r) {
    if (!r) throw PyErrorSet();
    return Expr(r);
}

struct Syms {
    Expr add, mul, integer, symbol, sympify, zero, one;
    Syms() {
        Expr mod = checked(PyImport_ImportModule("symengine"));
        add = checked(PyObject_GetAttrString(mod.get(), "Add"));
        mul = checked(PyObject_GetAttrString(mod.get(), "Mul"));
        integer = checked(PyObject_GetAttrString(mod.get(), "Integer"));
        symbol = checked(PyObject_GetAttrString(mod.get(), "Symbol"));
        sympify = checked(PyObject_GetAttrString(mod.get(), "sympify"));
        zero = checked(PyObject_CallFunction(integer.get(), "i", 0));
        one = checked(PyObject_CallFunction(integer.get(), "i", 1));
    }
};
// Heap-allocated and never freed: a static object's destructor would Py_DECREF after the
// interpreter has finalized.
inline Syms& syms() {
    static Syms* s = new Syms();
    return *s;
}

inline Expr call_star(const Expr& f, const ExprVec& args) {
    Expr tup = checked(PyTuple_New((Py_ssize_t)args.size()));
    for (size_t i = 0; i < args.size(); ++i) {
        Py_INCREF(args[i].get());
        PyTuple_SET_ITEM(tup.get(), (Py_ssize_t)i, args[i].get());
    }
    return checked(PyObject_CallObject(f.get(), tup.get()));
}
}  // namespace pyexpr_detail

inline Expr ex_zero() { return pyexpr_detail::syms().zero; }
inline Expr ex_one() { return pyexpr_detail::syms().one; }
inline Expr ex_integer(long n) { return pyexpr_detail::checked(PyObject_CallFunction(pyexpr_detail::syms().integer.get(), "l", n)); }
inline Expr ex_symbol(const std::string& name) { return pyexpr_detail::checked(PyObject_CallFunction(pyexpr_detail::syms().symbol.get(), "s", name.c_str())); }
inline Expr ex_add(const Expr& a, const Expr& b) { return pyexpr_detail::checked(PyNumber_Add(a.get(), b.get())); }
inline Expr ex_sub(const Expr& a, const Expr& b) { return pyexpr_detail::checked(PyNumber_Subtract(a.get(), b.get())); }
inline Expr ex_mul(const Expr& a, const Expr& b) { return pyexpr_detail::checked(PyNumber_Multiply(a.get(), b.get())); }
inline Expr ex_neg(const Expr& a) { return pyexpr_detail::checked(PyNumber_Negative(a.get())); }
inline Expr ex_pow(const Expr& a, int n) {
    Expr e = ex_integer(n);
    return pyexpr_detail::checked(PyNumber_Power(a.get(), e.get(), Py_None));
}
inline Expr ex_add(const ExprVec& v) {
    if (v.empty()) return ex_zero();
    if (v.size() == 1) return v[0];
    return pyexpr_detail::call_star(pyexpr_detail::syms().add, v);
}
inline Expr ex_mul(const ExprVec& v) {
    if (v.empty()) return ex_one();
    if (v.size() == 1) return v[0];
    return pyexpr_detail::call_star(pyexpr_detail::syms().mul, v);
}
inline bool ex_is_zero(const Expr& e) {
    int r = PyObject_RichCompareBool(e.get(), pyexpr_detail::syms().zero.get(), Py_EQ);
    if (r < 0) throw PyErrorSet();
    return r == 1;
}
inline Expr ex_sympify(const Expr& e) { return pyexpr_detail::checked(PyObject_CallFunctionObjArgs(pyexpr_detail::syms().sympify.get(), e.get(), nullptr)); }
inline std::string ex_str(const Expr& e) {
    Expr s = pyexpr_detail::checked(PyObject_Str(e.get()));
    const char* c = PyUnicode_AsUTF8(s.get());
    if (!c) throw PyErrorSet();
    return c;
}

#else
// ---------------------------------------------------------------------------
// SymEngine C++ backend
// ---------------------------------------------------------------------------
#include <symengine/add.h>
#include <symengine/basic.h>
#include <symengine/constants.h>
#include <symengine/integer.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/printers.h>
#include <symengine/symbol.h>

typedef SymEngine::RCP<const SymEngine::Basic> Expr;
typedef SymEngine::vec_basic ExprVec;

inline Expr ex_zero() { return SymEngine::zero; }
inline Expr ex_one() { return SymEngine::one; }
inline Expr ex_integer(long n) { return SymEngine::integer(n); }
inline Expr ex_symbol(const std::string& name) { return SymEngine::symbol(name); }
inline Expr ex_add(const Expr& a, const Expr& b) { return SymEngine::add(a, b); }
inline Expr ex_sub(const Expr& a, const Expr& b) { return SymEngine::sub(a, b); }
inline Expr ex_mul(const Expr& a, const Expr& b) { return SymEngine::mul(a, b); }
inline Expr ex_neg(const Expr& a) { return SymEngine::neg(a); }
inline Expr ex_pow(const Expr& a, int n) { return SymEngine::pow(a, SymEngine::integer(n)); }
inline Expr ex_add(const ExprVec& v) { return SymEngine::add(v); }
inline Expr ex_mul(const ExprVec& v) { return SymEngine::mul(v); }
inline bool ex_is_zero(const Expr& e) { return SymEngine::eq(*e, *SymEngine::zero); }
inline std::string ex_str(const Expr& e) { return e->__str__(); }
#endif
