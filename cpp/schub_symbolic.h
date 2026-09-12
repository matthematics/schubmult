// schub_symbolic.h: SymEngine layer shared by the double kernels — factorial elementary
// symmetric polynomials (schub_poly.elem_sym_poly), their memo, and call_zvars.

#pragma once

#include "schub_common.h"

#include <symengine/add.h>
#include <symengine/basic.h>
#include <symengine/constants.h>
#include <symengine/integer.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/printers.h>
#include <symengine/symbol.h>

using SymEngine::Basic;
using SymEngine::RCP;
using SymEngine::vec_basic;
typedef RCP<const Basic> Expr;

static bool is_zero(const Expr& e) { return SymEngine::eq(*e, *SymEngine::zero); }

// e_p(y[xs..xs+k) | z[ys..]) via the same divide-and-conquer recursion as the Python.
static Expr elem_sym_poly(int p, int k, const std::vector<Expr>& y, const std::vector<Expr>& z, int xs, int ys) {
    if (p > k) return SymEngine::zero;
    if (p == 0) return SymEngine::one;
    if (p == 1) {
        vec_basic terms;
        for (int i = 0; i < k; ++i) terms.push_back(SymEngine::sub(y.at(xs + i), z.at(ys + i)));
        return SymEngine::add(terms);
    }
    if (p == k) {
        Expr res = SymEngine::mul(SymEngine::sub(y.at(xs), z.at(ys)), SymEngine::sub(y.at(xs + 1), z.at(ys)));
        for (int i = 2; i < k; ++i) res = SymEngine::mul(res, SymEngine::sub(y.at(i + xs), z.at(ys)));
        return res;
    }
    int mid = k / 2, xsm = xs + mid, ysm = ys + mid, kmm = k - mid;
    Expr res = SymEngine::add(elem_sym_poly(p, mid, y, z, xs, ys), elem_sym_poly(p, kmm, y, z, xsm, ysm));
    for (int p2 = std::max(1, p - kmm); p2 < std::min(p, mid + 1); ++p2)
        res = SymEngine::add(res, SymEngine::mul(elem_sym_poly(p2, mid, y, z, xs, ys), elem_sym_poly(p - p2, kmm, y, z, xsm, ysm - p2)));
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
// y's and in the z's, so the key uses the sorted index lists.
struct ElemSymCache {
    std::vector<Expr> Y, Z, Q;  // Y[i] = y_i, Z[i] = z_i, Q[i] = q_i (1-indexed)
    std::unordered_map<std::vector<int>, Expr, VecHash> memo;
    size_t misses = 0;

    ElemSymCache(bool same) {
        Y.resize(MAXN + 1);
        Z.resize(MAXN + 1);
        Q.resize(MAXN + 1);
        for (int i = 0; i <= MAXN; ++i) {
            Y[i] = SymEngine::symbol("y_" + std::to_string(i));
            Z[i] = same ? Y[i] : Expr(SymEngine::symbol("z_" + std::to_string(i)));
            Q[i] = SymEngine::symbol("q_" + std::to_string(i));
        }
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
        for (int i : yidx) yv.push_back(Y.at(i));
        for (int i : zidx) zv.push_back(Z.at(i));
        Expr e = elem_sym_poly(p, k, yv, zv, 0, 0);
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
