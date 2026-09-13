// kernel_double_alt.h: C++ port of schubmult.mult.double.schubmult_double_alt_from_elems_backwards
// (= schubmult_double_alt_from_elems), the pull-out-variable recursion used for positivity,
// together with schub_lib.pull_out_var and elem_sym_positional_perms.
//
// The elementary symmetric objects come from esc.builder, so this reproduces the Python
// exactly for any elem_func (FactorialElemSym, elem_sym_poly, ...).

#pragma once

#include "schub_symbolic.h"

#include <map>

// Python's len() of a trimmed Permutation: max(trimmed length, 2). Only used where the
// Python algorithm's loop bounds depend on it; the representation stays fixed-size.
static int py_len(const Perm& a) {
    int L = MAXN;
    while (L > 0 && a.p[L - 1] == L) --L;
    return std::max(L, 2);
}

static int perm_inv_full(const Perm& a) { return perm_inv(a, MAXN); }

struct PullOut {
    std::vector<int> index_list;  // 1-based values of v at the retained positions
    Perm vp;
    bool operator<(const PullOut& o) const {
        if (vp != o.vp) return vp < o.vp;
        return index_list < o.index_list;
    }
    bool operator==(const PullOut& o) const { return vp == o.vp && index_list == o.index_list; }
};

// schub_lib.pull_out_var(vnum, v): the set of (index_list, vp).
static std::vector<PullOut> pull_out_var(int vnum, const Perm& v) {
    int L = py_len(v);
    std::vector<PullOut> ret;
    if (vnum >= L) {
        ret.push_back({{}, v});
        return ret;
    }
    if (L + 2 > MAXN) die("pull_out_var: permutation too long for MAXN");
    auto add_result = [&](const Perm& vpm) {
        // vpm with position vnum-1 (holding L+1) removed, values above L+1 shifted down
        Perm vp = identity_perm();
        int out = 0;
        int lenv = py_len(vpm);
        for (int i = 0; i < lenv; ++i) {
            if (i == vnum - 1) continue;
            int x = vpm.p[i];
            vp.p[out++] = (uint8_t)(x > L + 1 ? x - 1 : x);
        }
        int lvp = py_len(vp);
        PullOut po;
        po.vp = vp;
        for (int i = vnum; i < L; ++i)
            if ((i > lvp && v.p[i] == i) || (i <= lvp && v.p[i] == vp.p[i - 1])) po.index_list.push_back(v.p[i]);
        ret.push_back(std::move(po));
    };
    std::vector<std::pair<Perm, int>> cur, nxt;
    cur.push_back({v, 0});
    for (int p = 0; p < L + 1 - vnum; ++p) {
        nxt.clear();
        for (const auto& e : cur) {
            const Perm& vpm = e.first;
            int b = e.second;
            if (vpm.p[vnum - 1] == L + 1) add_result(vpm);
            for (int j = vnum; j < L + 2; ++j) {
                if (vpm.p[j] <= b) continue;
                for (int i = 0; i < vnum; ++i)
                    if (bruhat_ascent(vpm, i, j)) nxt.push_back({swapped(vpm, i, j), vpm.p[j]});
            }
        }
        std::swap(cur, nxt);
    }
    for (const auto& e : cur)
        if (e.first.p[vnum - 1] == L + 1) add_result(e.first);
    std::sort(ret.begin(), ret.end());
    ret.erase(std::unique(ret.begin(), ret.end()), ret.end());
    return ret;
}

struct PosUp {
    Perm perm;
    int udiff;
    int sign;
    bool operator<(const PosUp& o) const {
        if (perm != o.perm) return perm < o.perm;
        if (udiff != o.udiff) return udiff < o.udiff;
        return sign < o.sign;
    }
    bool operator==(const PosUp& o) const { return perm == o.perm && udiff == o.udiff && sign == o.sign; }
};

// schub_lib.elem_sym_positional_perms(orig, p, *k) with k given 0-based; returns the set of
// (perm, udiff, sign).
static std::vector<PosUp> elem_sym_positional_perms(const Perm& orig, int p, const std::vector<int>& kpos) {
    std::vector<PosUp> total{{orig, 0, 1}};
    std::vector<std::pair<Perm, int>> cur{{orig, 1}}, nxt;
    int maxk = 0;
    for (int q : kpos) maxk = std::max(maxk, q);
    bool ink[MAXN] = {false};
    for (int q : kpos) ink[q] = true;
    for (int pp = 0; pp < p; ++pp) {
        nxt.clear();
        for (const auto& e : cur) {
            const Perm& up = e.first;
            int sign = e.second;
            int bound = py_len(up) + maxk + 1;
            if (bound > MAXN) die("elem_sym_positional_perms: permutation too long for MAXN");
            for (int i : kpos) {
                if (up.p[i] != orig.p[i]) continue;
                for (int j = 0; j < bound; ++j) {
                    if (ink[j]) continue;
                    int a = std::min(i, j), b = std::max(i, j);
                    if (!bruhat_ascent(up, a, b)) continue;
                    Perm np = swapped(up, a, b);
                    int ns = i < j ? sign : -sign;
                    nxt.push_back({np, ns});
                    total.push_back({np, pp + 1, ns});
                }
            }
        }
        std::sort(nxt.begin(), nxt.end());
        nxt.erase(std::unique(nxt.begin(), nxt.end()), nxt.end());
        std::swap(cur, nxt);
    }
    std::sort(total.begin(), total.end());
    total.erase(std::unique(total.begin(), total.end()), total.end());
    return total;
}

#ifndef SCHUB_EXPRDICT
#define SCHUB_EXPRDICT
typedef std::vector<std::pair<Perm, Expr>> ExprDict;
#endif

typedef std::unordered_map<Perm, vec_basic, PermHash> TermMap;

static TermMap alt_backwards_rec(const TermMap& perm_dict, const Perm& v, ElemSymCache& esc) {
    if (perm_inv_full(v) == 0) return perm_dict;
    Perm vinv = inverse(v, MAXN);
    int index = 0;  // max descent of ~v, 1-based
    for (int i = 0; i + 1 < MAXN; ++i)
        if (vinv.p[i] > vinv.p[i + 1]) index = i + 1;
    std::vector<PullOut> L = pull_out_var(index, vinv);

    TermMap ret;
    std::map<Perm, TermMap> cache;  // by new_v, as the Python
    std::vector<int> kpos, yidx;
    for (const PullOut& po : L) {
        auto it = cache.find(po.vp);
        if (it == cache.end()) it = cache.emplace(po.vp, alt_backwards_rec(perm_dict, inverse(po.vp, MAXN), esc)).first;
        const TermMap& start = it->second;
        int m = (int)po.index_list.size();
        kpos.clear();
        for (int i : po.index_list) kpos.push_back(i - 1);
        for (const auto& kv : start) {
            const Perm& u = kv.first;
            Expr val = SymEngine::add(kv.second);
            if (is_zero(val)) continue;
            for (const PosUp& e : elem_sym_positional_perms(u, m, kpos)) {
                yidx.clear();
                for (int i : po.index_list)
                    if (e.perm.p[i - 1] == u.p[i - 1]) yidx.push_back(e.perm.p[i - 1]);
                const Expr& ef = esc.get(m - e.udiff, m - e.udiff, yidx, std::vector<int>{index});
                if (is_zero(ef)) continue;
                Expr term = SymEngine::mul(val, ef);
                if (e.sign < 0) term = SymEngine::neg(term);
                ret[e.perm].push_back(term);
            }
        }
    }
    return ret;
}

// schubmult_double_alt_from_elems_backwards(perm_dict, v, var2, var3, elem_func)
static ExprDict schubmult_double_alt_from_elems_backwards(const ExprDict& perm_dict, const Perm& v, ElemSymCache& esc) {
    TermMap in;
    for (const auto& kv : perm_dict) in[kv.first].push_back(kv.second);
    TermMap out = alt_backwards_rec(in, v, esc);
    ExprDict r;
    for (auto& kv : out) {
        Expr s = SymEngine::add(kv.second);
        if (!is_zero(s)) r.push_back({kv.first, s});
    }
    std::sort(r.begin(), r.end(), [](const std::pair<Perm, Expr>& x, const std::pair<Perm, Expr>& y) { return x.first < y.first; });
    return r;
}
