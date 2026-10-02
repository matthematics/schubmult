// kernel_double_transition.h: double Schubert products by the double transition formula + double Monk.
//
// Transition for double Schubert polynomials (Kohnert--Veigneau, Prop. 4.1; Lascoux): with r the last
// descent of w, s the last position with w(s) < w(r) and v = w t_{rs},
//     S_w(x, z) = (x_r - z_{v(r)}) S_v(x, z) + sum_{i < r, v <. v t_{ir}} S_{v t_{ir}}(x, z).
// For a product (sum_u c_u S_u(x, y)) * S_v(x, z) the recursion is run on the dict itself, P(w) =
// S_w(x, z) * D: the linear factor acts by the double Monk rule
//     (x_r - z_c) S_u(x, y) = sum_{j > r} S_{u t_{rj}} - sum_{i < r} S_{u t_{ir}} + (y_{u(r)} - z_c) S_u
// (covers of u only), and the recursion is memoized on w, which is the sharing lrcalc's Horner scheme
// provides in the single case.  Coefficients stay structural (products of y - z binomials summed per
// permutation); nothing is expanded.

#pragma once

#include <memory>

#include "schub_symbolic.h"

#ifndef SCHUB_EXPRDICT
#define SCHUB_EXPRDICT
typedef std::vector<std::pair<Perm, Expr>> ExprDict;
#endif

namespace dtrans_detail {

typedef std::unordered_map<Perm, Expr, PermHash> ExprMap;
typedef std::unordered_map<Perm, ExprVec, PermHash> PendingMap;

static void pend(PendingMap& acc, const Perm& w, Expr e) { acc[w].push_back(std::move(e)); }

// acc += (x_r - z_c) * sum_u in[u] S_u(x, y), in S_n.
static void linear_factor(int r, const Expr& zc, const ExprMap& in, int n, ElemSymCache& esc, PendingMap& acc) {
    for (const auto& kv : in) {
        const Perm& u = kv.first;
        const Expr& c = kv.second;
        int ur = u.p[r - 1];
        int last = 0;
        for (int j = r - 1; j >= 1; --j) {
            int uj = u.p[j - 1];
            if (last < uj && uj < ur) {
                last = uj;
                Perm t = u;
                std::swap(t.p[j - 1], t.p[r - 1]);
                pend(acc, t, ex_neg(c));
            }
        }
        last = n + 1;
        for (int j = r + 1; j <= n; ++j) {
            int uj = u.p[j - 1];
            if (ur < uj && uj < last) {
                last = uj;
                Perm t = u;
                std::swap(t.p[r - 1], t.p[j - 1]);
                pend(acc, t, c);
            }
        }
        Expr f = ex_sub(ElemSymCache::need(esc.Y, ur), zc);
        if (!is_zero(f)) pend(acc, u, ex_mul(c, f));
    }
}

static std::shared_ptr<const ExprMap> collapse(PendingMap& acc) {
    auto out = std::make_shared<ExprMap>();
    for (auto& kv : acc) {
        Expr s = ex_add(kv.second);
        if (!is_zero(s)) out->emplace(kv.first, std::move(s));
    }
    return out;
}

struct Rec {
    const std::shared_ptr<const ExprMap> base;
    int n;
    ElemSymCache& esc;
    std::unordered_map<Perm, std::shared_ptr<const ExprMap>, PermHash> memo;

    std::shared_ptr<const ExprMap> run(const Perm& w) {
        int r = 0;
        for (int i = n - 1; i >= 1; --i)
            if (w.p[i - 1] > w.p[i]) { r = i; break; }
        if (r == 0) return base;
        auto it = memo.find(w);
        if (it != memo.end()) return it->second;
        int s = r + 1;
        while (s < n && w.p[r - 1] > w.p[s]) ++s;
        Perm v = w;
        std::swap(v.p[r - 1], v.p[s - 1]);
        PendingMap acc;
        linear_factor(r, ElemSymCache::need(esc.Z, v.p[r - 1]), *run(v), n, esc, acc);
        int vr = v.p[r - 1], last = 0;
        for (int i = r - 1; i >= 1; --i) {
            int vi = v.p[i - 1];
            if (last < vi && vi < vr) {
                last = vi;
                Perm next = v;
                std::swap(next.p[i - 1], next.p[r - 1]);
                for (const auto& kv : *run(next)) pend(acc, kv.first, kv.second);
            }
        }
        return memo.emplace(w, collapse(acc)).first->second;
    }
};

}  // namespace dtrans_detail

// (sum_u c_u S_u(x, y)) * S_v(x, z) in the basis S_w(x, y), by the double transition recursion on v.
static ExprDict schubmult_double_transition(const ExprDict& perm_dict, const Perm& v, int n, ElemSymCache& esc) {
    auto base = std::make_shared<dtrans_detail::ExprMap>();
    for (const auto& kv : perm_dict)
        if (!is_zero(kv.second)) (*base)[kv.first] = kv.second;
    dtrans_detail::Rec rec{base, n, esc, {}};
    std::shared_ptr<const dtrans_detail::ExprMap> res = rec.run(v);
    ExprDict out(res->begin(), res->end());
    std::sort(out.begin(), out.end(), [](const std::pair<Perm, Expr>& x, const std::pair<Perm, Expr>& y) { return x.first < y.first; });
    return out;
}
