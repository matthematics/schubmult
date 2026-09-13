// schubmult_api.h: C++ entry points for the schubmult_cpp extension (schubmult_module.cpp).
// Compiled with SCHUB_PYEXPR: coefficients are Python `symengine` objects (see expr.h), so this
// header has no SymEngine C++ dependency. Permutations cross the boundary as one-line-notation
// int vectors and generating sets as vectors of Python objects (index i = i-th variable).

#pragma once

#define SCHUB_THROW
#include "kernel_double.h"
#include "kernel_double_alt.h"
#include "kernel_q.h"
#include "kernel_q_double.h"
#include "kernel_single.h"

#include <map>
#include <memory>

typedef std::vector<int> PyPerm;
typedef std::vector<std::pair<PyPerm, int>> PyIntDict;
typedef std::vector<std::pair<PyPerm, Expr>> PyExprDict;
typedef std::vector<Expr> PyVars;

static Perm perm_from_py(const PyPerm& arr) {
    int L = (int)arr.size();
    if (L > MAXN) die("permutation too long for MAXN; rebuild the extension with a larger MAXN");
    std::vector<char> seen(L + 1, 0);
    for (int x : arr) {
        if (x < 1 || x > L || seen[x]) die("input is not a permutation of 1..n");
        seen[x] = 1;
    }
    Perm r = identity_perm();
    for (int i = 0; i < L; ++i) r.p[i] = (uint8_t)arr[i];
    return r;
}

static PyPerm perm_to_py(const Perm& a) {
    int L = MAXN;
    while (L > 0 && a.p[L - 1] == L) --L;
    return PyPerm(a.p, a.p + L);
}

static int perm_len(const PyPerm& a) {
    int L = (int)a.size();
    while (L > 0 && a[L - 1] == L) --L;
    return L;
}

static std::string vars_key(const PyVars& v) {
    std::string s;
    for (const Expr& e : v) {
        s += e.is_null() ? std::string("~") : ex_str(e);
        s += '\x01';
    }
    return s;
}

// Builds e_p(y | z) by calling a Python elem_func(p, k, yvars, zvars) (the *_from_elems kernels).
static ElemSymCache::Builder python_builder(const Expr& elem_func) {
    return [elem_func](int p, int k, const std::vector<Expr>& y, const std::vector<Expr>& z) {
        auto list_of = [](const std::vector<Expr>& v) {
            Expr lst = pyexpr_detail::checked(PyList_New((Py_ssize_t)v.size()));
            for (size_t i = 0; i < v.size(); ++i) {
                Py_INCREF(v[i].get());
                PyList_SET_ITEM(lst.get(), (Py_ssize_t)i, v[i].get());
            }
            return lst;
        };
        Expr yl = list_of(y), zl = list_of(z);
        Expr r = pyexpr_detail::checked(PyObject_CallFunction(elem_func.get(), "iiOO", p, k, yl.get(), zl.get()));
        return ex_sympify(r);
    };
}

// One ElemSymCache per distinct (y, z, q, elem_func), so the e_p memo survives across calls.
static ElemSymCache& esc_for(const PyVars& y, const PyVars& z, const PyVars& q, const Expr& elem_func = Expr()) {
    // leaked on purpose: the caches hold Python references (see expr.h)
    static auto* caches = new std::map<std::string, std::unique_ptr<ElemSymCache>>();
    std::string key = vars_key(y) + "\x02" + vars_key(z) + "\x02" + vars_key(q) + "\x02" + std::to_string((uintptr_t)elem_func.get());
    auto it = caches->find(key);
    if (it == caches->end()) {
        it = caches->emplace(key, std::make_unique<ElemSymCache>(y, z, q)).first;
        if (!elem_func.is_null()) it->second->builder = python_builder(elem_func);  // keeps elem_func alive with the cache
    }
    return *it->second;
}

static int dict_bound(const PyExprDict& d, const PyPerm& vpy) {
    int a_max = 1;
    for (const auto& kv : d) a_max = std::max(a_max, perm_len(kv.first));
    int n = std::max(2, a_max + std::max(1, perm_len(vpy)) - 1);  // S_a * S_b lives in S_{a+b-1}
    if (n > MAXN) die("product needs S_n with n > MAXN; rebuild the extension with a larger MAXN");
    return n;
}

static ExprDict expr_dict_in(const PyExprDict& d) {
    ExprDict in;
    for (const auto& kv : d) in.push_back({perm_from_py(kv.first), kv.second});
    return in;
}

static PyExprDict expr_dict_out(const ExprDict& out) {
    PyExprDict r;
    for (const auto& kv : out) r.push_back({perm_to_py(kv.first), kv.second});
    return r;
}

static PyIntDict api_schubmult_py(const PyIntDict& d, const PyPerm& vpy) {
    IntDict in;
    int a_max = 1;
    for (const auto& kv : d) {
        in.push_back({perm_from_py(kv.first), kv.second});
        a_max = std::max(a_max, perm_len(kv.first));
    }
    int n = std::max(2, a_max + std::max(1, perm_len(vpy)) - 1);
    if (n > MAXN) die("product needs S_n with n > MAXN; rebuild the extension with a larger MAXN");
    IntDict out = schubmult_single(in, perm_from_py(vpy), n);
    PyIntDict r;
    for (const auto& kv : out) r.push_back({perm_to_py(kv.first), kv.second});
    return r;
}

static PyExprDict api_schubmult_double(const PyExprDict& d, const PyPerm& vpy, const PyVars& y, const PyVars& z, const Expr& elem_func) {
    return expr_dict_out(schubmult_double(expr_dict_in(d), perm_from_py(vpy), dict_bound(d, vpy), esc_for(y, z, PyVars(), elem_func)));
}

static PyExprDict api_schubmult_double_alt_from_elems(const PyExprDict& d, const PyPerm& vpy, const PyVars& y, const PyVars& z, const Expr& elem_func) {
    return expr_dict_out(schubmult_double_alt_from_elems_backwards(expr_dict_in(d), perm_from_py(vpy), esc_for(y, z, PyVars(), elem_func)));
}

// Integer input coefficients; output is symbolic in the q's.
static PyExprDict api_schubmult_q(const PyIntDict& d, const PyPerm& vpy, const PyVars& q) {
    ElemSymCache& esc = esc_for(PyVars(), PyVars(), q);
    QDict in;
    for (const auto& kv : d) {
        QPoly c;
        c.t.push_back({mono_one(), kv.second});
        in.push_back({perm_from_py(kv.first), c});
    }
    QDict out = schubmult_q_fast(in, perm_from_py(vpy));
    PyExprDict r;
    for (const auto& kv : out) {
        ExprVec terms;
        for (const auto& t : kv.second.t) terms.push_back(ex_mul(ex_integer(t.second), mono_expr(t.first, esc)));
        r.push_back({perm_to_py(kv.first), ex_add(terms)});
    }
    return r;
}

static PyExprDict api_schubmult_q_double(const PyExprDict& d, const PyPerm& vpy, const PyVars& y, const PyVars& z, const PyVars& q) {
    return expr_dict_out(schubmult_q_double_fast(expr_dict_in(d), perm_from_py(vpy), esc_for(y, z, q)));
}
