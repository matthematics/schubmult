// schubmult_api.h: thin C++ entry points for the Python extension (schubmult_cpp.pyx).
// Permutations cross the boundary as one-line-notation int vectors, coefficients as
// SymEngine RCP<const Basic> (or int for the classical kernel), and generating sets as
// vectors of expressions (index i = i-th variable), so any symbols the caller uses work.

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
        s += e.is_null() ? std::string("~") : e->__str__();
        s += '\x01';
    }
    return s;
}

// One ElemSymCache per distinct (y, z, q, elem_func) so the e_p memo survives across calls.
// `cb(ctx, p, k, y, z)` builds the elementary symmetric object (null ctx: elem_sym_poly).
typedef Expr (*ElemFuncCB)(void* ctx, int p, int k, const std::vector<Expr>& y, const std::vector<Expr>& z);

static ElemSymCache& esc_for(const PyVars& y, const PyVars& z, const PyVars& q, ElemFuncCB cb = nullptr, void* ctx = nullptr) {
    static std::map<std::string, std::unique_ptr<ElemSymCache>> caches;
    std::string key = vars_key(y) + "\x02" + vars_key(z) + "\x02" + vars_key(q) + "\x02" + std::to_string((uintptr_t)ctx);
    auto it = caches.find(key);
    if (it == caches.end()) {
        it = caches.emplace(key, std::make_unique<ElemSymCache>(y, z, q)).first;
        if (ctx) it->second->builder = [cb, ctx](int p, int k, const std::vector<Expr>& yv, const std::vector<Expr>& zv) {
            Expr e = cb(ctx, p, k, yv, zv);
            if (e.is_null()) die("elem_func raised");
            return e;
        };
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

static PyIntDict api_schubmult_py(const PyIntDict& d, const PyPerm& vpy) {
    IntDict in;
    int a_max = 1;
    for (const auto& kv : d) {
        in.push_back({perm_from_py(kv.first), kv.second});
        a_max = std::max(a_max, perm_len(kv.first));
    }
    int n = std::max(2, a_max + std::max(1, perm_len(vpy)) - 1);  // S_a * S_b lives in S_{a+b-1}
    if (n > MAXN) die("product needs S_n with n > MAXN; rebuild the extension with a larger MAXN");
    IntDict out = schubmult_single(in, perm_from_py(vpy), n);
    PyIntDict r;
    for (const auto& kv : out) r.push_back({perm_to_py(kv.first), kv.second});
    return r;
}

static PyExprDict api_schubmult_double(const PyExprDict& d, const PyPerm& vpy, const PyVars& y, const PyVars& z) {
    ExprDict in;
    for (const auto& kv : d) in.push_back({perm_from_py(kv.first), kv.second});
    ExprDict out = schubmult_double(in, perm_from_py(vpy), dict_bound(d, vpy), esc_for(y, z, PyVars()));
    PyExprDict r;
    for (const auto& kv : out) r.push_back({perm_to_py(kv.first), kv.second});
    return r;
}

// schubmult_double_from_elems: the theta-code kernel with elem_func in place of elem_sym_poly.
static PyExprDict api_schubmult_double_from_elems(const PyExprDict& d, const PyPerm& vpy, const PyVars& y, const PyVars& z, ElemFuncCB cb, void* ctx) {
    ExprDict in;
    for (const auto& kv : d) in.push_back({perm_from_py(kv.first), kv.second});
    ExprDict out = schubmult_double(in, perm_from_py(vpy), dict_bound(d, vpy), esc_for(y, z, PyVars(), cb, ctx));
    PyExprDict r;
    for (const auto& kv : out) r.push_back({perm_to_py(kv.first), kv.second});
    return r;
}

static PyExprDict api_schubmult_double_alt_from_elems(const PyExprDict& d, const PyPerm& vpy, const PyVars& y, const PyVars& z, ElemFuncCB cb, void* ctx) {
    ExprDict in;
    for (const auto& kv : d) in.push_back({perm_from_py(kv.first), kv.second});
    ExprDict out = schubmult_double_alt_from_elems_backwards(in, perm_from_py(vpy), esc_for(y, z, PyVars(), cb, ctx));
    PyExprDict r;
    for (const auto& kv : out) r.push_back({perm_to_py(kv.first), kv.second});
    return r;
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
        vec_basic terms;
        for (const auto& t : kv.second.t) terms.push_back(SymEngine::mul(SymEngine::integer(t.second), mono_expr(t.first, esc)));
        r.push_back({perm_to_py(kv.first), SymEngine::add(terms)});
    }
    return r;
}

static PyExprDict api_schubmult_q_double(const PyExprDict& d, const PyPerm& vpy, const PyVars& y, const PyVars& z, const PyVars& q) {
    ExprDict in;
    for (const auto& kv : d) in.push_back({perm_from_py(kv.first), kv.second});
    ExprDict out = schubmult_q_double_fast(in, perm_from_py(vpy), esc_for(y, z, q));
    PyExprDict r;
    for (const auto& kv : out) r.push_back({perm_to_py(kv.first), kv.second});
    return r;
}
