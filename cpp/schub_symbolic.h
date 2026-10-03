// schub_symbolic.h: symbolic layer shared by the double kernels (see expr.h for the backends) — factorial elementary
// symmetric polynomials (schub_poly.elem_sym_poly), their memo, and call_zvars.

#pragma once

#include "schub_common.h"

#include <functional>
#include <unordered_map>

#include "expr.h"

#ifdef _MSC_VER
#include <intrin.h>
#endif

static bool is_zero(const Expr& e) { return ex_is_zero(e); }

// ---------------------------------------------------------------------------
// Shadow values: a partial sum's value at SHADOW_POINTS random points of F_p, p = 2^61 - 1,
// carried alongside the expression and combined node by node (see schubmult/mult/_shadow.py).
// A state whose shadow vanishes at every point is dropped (probabilistic=True).
// ---------------------------------------------------------------------------

static const uint64_t SHADOW_P = (1ULL << 61) - 1;
static const int SHADOW_POINTS = 2;

static inline uint64_t shadow_mulmod(uint64_t a, uint64_t b) {
#ifdef _MSC_VER
    uint64_t hi, lo = _umul128(a, b, &hi);
#else
    unsigned __int128 m = (unsigned __int128)a * b;
    uint64_t lo = (uint64_t)m, hi = (uint64_t)(m >> 64);
#endif
    // 2^61 = 1 mod p, so a*b = hi*2^64 + lo = 8*hi + (lo >> 61) + (lo & p)
    uint64_t r = (lo & SHADOW_P) + ((lo >> 61) | (hi << 3));
    r = (r & SHADOW_P) + (r >> 61);
    return r >= SHADOW_P ? r - SHADOW_P : r;
}

struct ShadowVal {
    uint64_t v[SHADOW_POINTS] = {0, 0};
    bool is_zero() const {
        for (int i = 0; i < SHADOW_POINTS; ++i)
            if (v[i]) return false;
        return true;
    }
    ShadowVal operator*(const ShadowVal& o) const {
        ShadowVal r;
        for (int i = 0; i < SHADOW_POINTS; ++i) r.v[i] = shadow_mulmod(v[i], o.v[i]);
        return r;
    }
    ShadowVal operator-() const {
        ShadowVal r;
        for (int i = 0; i < SHADOW_POINTS; ++i) r.v[i] = v[i] ? SHADOW_P - v[i] : 0;
        return r;
    }
    ShadowVal& operator+=(const ShadowVal& o) {
        for (int i = 0; i < SHADOW_POINTS; ++i) {
            uint64_t s = v[i] + o.v[i];
            v[i] = s >= SHADOW_P ? s - SHADOW_P : s;
        }
        return *this;
    }
};

// Evaluates kernel *inputs* (coefficients, e_p's) at the sample points; memoized by object
// identity, which is stable for the lifetime of a call (inputs live in the caller's dict and in
// the ElemSymCache). Null `eval` means exact mode.
struct Shadow {
    std::function<ShadowVal(const Expr&)> eval;
    std::unordered_map<const void*, ShadowVal> memo;

    const ShadowVal& of(const Expr& e) {
        auto it = memo.find((const void*)e.get());
        if (it != memo.end()) return it->second;
        return memo.emplace((const void*)e.get(), eval(e)).first->second;
    }
};

// e_p(y[xs..xs+k) | z[ys..]) via the same divide-and-conquer recursion as the Python.
static Expr elem_sym_poly(int p, int k, const std::vector<Expr>& y, const std::vector<Expr>& z, int xs, int ys) {
    if (p > k) return ex_zero();
    if (p == 0) return ex_one();
    if (p == 1) {
        ExprVec terms;
        for (int i = 0; i < k; ++i) terms.push_back(ex_sub(y.at(xs + i), z.at(ys + i)));
        return ex_add(terms);
    }
    if (p == k) {
        Expr res = ex_mul(ex_sub(y.at(xs), z.at(ys)), ex_sub(y.at(xs + 1), z.at(ys)));
        for (int i = 2; i < k; ++i) res = ex_mul(res, ex_sub(y.at(i + xs), z.at(ys)));
        return res;
    }
    int mid = k / 2, xsm = xs + mid, ysm = ys + mid, kmm = k - mid;
    Expr res = ex_add(elem_sym_poly(p, mid, y, z, xs, ys), elem_sym_poly(p, kmm, y, z, xsm, ysm));
    for (int p2 = std::max(1, p - kmm); p2 < std::min(p, mid + 1); ++p2)
        res = ex_add(res, ex_mul(elem_sym_poly(p2, mid, y, z, xs, ys), elem_sym_poly(p - p2, kmm, y, z, xsm, ysm - p2)));
    return res;
}

struct VecHash {
    size_t operator()(const std::vector<int>& v) const {
        uint64_t h = 0x9E3779B97F4A7C15ULL;
        for (int x : v) {
            h ^= (uint64_t)(uint32_t)x;
            h *= 0xFF51AFD7ED558CCDULL;
            h ^= h >> 33;
        }
        return (size_t)h;
    }
};

// Memoized e_p(y_{yidx} | z_{zidx}), kept unexpanded. The polynomial is symmetric in the
// y's and in the z's, so the key uses the sorted index lists. `builder` produces the
// elementary symmetric object (default: schub_poly.elem_sym_poly; the Python binding can
// substitute e.g. FactorialElemSym, the `elem_func` of the *_from_elems kernels).
struct ElemSymCache {
    typedef std::function<Expr(int p, int k, const std::vector<Expr>& y, const std::vector<Expr>& z)> Builder;
    std::vector<Expr> Y, Z, Q;  // Y[i] = y_i, Z[i] = z_i, Q[i] = q_i (1-indexed)
    std::unordered_map<std::vector<int>, Expr, VecHash> memo;
    Builder builder;
    size_t misses = 0;

    // yname/zname: symbol name prefixes ("y" -> y_1, y_2, ...); same variables when zname is empty.
    ElemSymCache(bool same, const std::string& yname = "y", const std::string& zname = "z", const std::string& qname = "q") {
        Y.resize(MAXN + 1);
        Z.resize(MAXN + 1);
        Q.resize(MAXN + 1);
        for (int i = 0; i <= MAXN; ++i) {
            Y[i] = ex_symbol(yname + "_" + std::to_string(i));
            Z[i] = same ? Y[i] : Expr(ex_symbol(zname + "_" + std::to_string(i)));
            Q[i] = ex_symbol(qname + "_" + std::to_string(i));
        }
    }

    // Explicit variable vectors (index i = i-th variable, 0..MAXN); any expression is allowed,
    // e.g. the integer 0 for a "zero" generating set. Missing entries (null) die when used.
    ElemSymCache(std::vector<Expr> y, std::vector<Expr> z, std::vector<Expr> q) : Y(std::move(y)), Z(std::move(z)), Q(std::move(q)) {
        Y.resize(MAXN + 1);
        Z.resize(MAXN + 1);
        Q.resize(MAXN + 1);
    }

    static const Expr& need(const std::vector<Expr>& v, int i) {
        if (v[i].is_null()) die("generating set has too few variables for this product");
        return v[i];
    }

    const Expr& get(int p, int k, std::vector<int> yidx, std::vector<int> zidx) {
        std::sort(yidx.begin(), yidx.end());
        std::sort(zidx.begin(), zidx.end());
        std::vector<int> key;
        key.reserve(yidx.size() + zidx.size() + 3);
        key.push_back(p);
        key.push_back(k);
        key.insert(key.end(), yidx.begin(), yidx.end());
        key.push_back(-1);
        key.insert(key.end(), zidx.begin(), zidx.end());
        auto it = memo.find(key);
        if (it != memo.end()) return it->second;
        ++misses;
        std::vector<Expr> yv, zv;
        for (int i : yidx) yv.push_back(need(Y, i));
        for (int i : zidx) zv.push_back(need(Z, i));
        Expr e = builder ? builder(p, k, yv, zv) : elem_sym_poly(p, k, yv, zv, 0, 0);
        return memo.emplace(key, e).first->second;
    }
};

// z-variable indices of elem_sym_func for the transition v1 -> v2 at 1-based level i
// (call_zvars): v2[i-1] followed by v2[j] over the positions j != i-1 where v1, v2 differ.
static std::vector<int> call_zvars(const Perm& v1, const Perm& v2, int i) {
    std::vector<int> r;
    r.push_back(v2.p[i - 1]);
    for (int j = 0; j < MAXN; ++j)
        if (j != i - 1 && v1.p[j] != v2.p[j]) r.push_back(v2.p[j]);
    return r;
}

// zidx[index][vpid][t]: the vdiff+1 z indices used by vp.trans[index][vpid][t]
// (elem_sym_poly(newk - vdiff, newk, ...) reads exactly that many).
typedef std::vector<std::vector<std::vector<std::vector<int>>>> ZIdx;
static ZIdx compute_zidx(const VPaths& vp) {
    int thL = (int)vp.trans.size();
    ZIdx zidx(thL);
    for (int index = 0; index < thL; ++index) {
        zidx[index].resize(vp.trans[index].size());
        for (size_t vpid = 0; vpid < vp.trans[index].size(); ++vpid)
            for (const Trans& tr : vp.trans[index][vpid]) {
                std::vector<int> zs = call_zvars(vp.level[index][vpid], vp.level[index + 1][tr.v2], index + 1);
                if ((int)zs.size() < tr.vdiff + 1) die("internal: too few z variables for a v-path step");
                zs.resize(tr.vdiff + 1);
                zidx[index][vpid].push_back(std::move(zs));
            }
    }
    return zidx;
}

// y indices of elem_sym_func(_q): y_{u2[j]} at the positions j < k where u1 and u2 agree.
static void yvars_of(const Perm& u1, const Perm& u2, int k, std::vector<int>& yidx) {
    yidx.clear();
    for (int j = 0; j < k; ++j)
        if (u1.p[j] == u2.p[j]) yidx.push_back(u2.p[j]);
}
