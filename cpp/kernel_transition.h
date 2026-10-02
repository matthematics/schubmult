// kernel_transition.h: single Schubert products by Lascoux--Schuetzenberger transition + Monk.
//
// This is the algorithm of Buch's lrcalc (schublib.c): expand one factor into monomials with the
// transition recursion S_w = x_r S_v + sum_i S_{v t_{ir}}, then multiply the monomials into the
// other factor one variable at a time by Monk's rule (Horner scheme over the Schubert basis).
// Its cost is governed by the number of pipe dreams of the expanded factor and is independent of
// the v-path structure, so it complements kernel_single.h: cheap exactly where the v-path kernel's
// signed intermediate expansion is large (many descents, large theta entries) and vice versa.
// schubmult_single_hybrid dispatches between the two by running the transition expansion with a
// pipe-dream budget.

#pragma once

#include <climits>
#include <cmath>
#include <map>

#include "kernel_single.h"

typedef std::vector<int> ExpMono;  // exponent vector, trailing zeros trimmed
typedef std::map<ExpMono, int_coef_t> ExpMonoDict;

struct TransBudgetExceeded {};

namespace trans_detail {

static void trim(ExpMono& m) {
    while (!m.empty() && m.back() == 0) m.pop_back();
}

static void add_mono(ExpMonoDict& out, ExpMono m, int_coef_t c) {
    trim(m);
    auto it = out.find(m);
    if (it == out.end()) {
        if (c) out.emplace(std::move(m), c);
    } else if ((it->second += c) == 0) {
        out.erase(it);
    }
}

// Transition recursion on w restricted to S_n; every leaf is one pipe dream, so `leaves` counts
// S_w(1,...,1) and lets a caller abort once the expansion is known to be large.
static void trans_rec(const Perm& w, int n, ExpMonoDict& out, long& leaves, long budget) {
    int r = 0;  // last descent, 1-indexed
    for (int i = n - 1; i >= 1; --i)
        if (w.p[i - 1] > w.p[i]) { r = i; break; }
    if (r == 0) {
        if (++leaves > budget) throw TransBudgetExceeded{};
        add_mono(out, ExpMono(), 1);
        return;
    }
    // s: the position after r carrying the largest value below w(r); positions after the last
    // descent are increasing, so scan while w(s) < w(r)
    int s = r + 1;
    while (s < n && w.p[r - 1] > w.p[s]) ++s;
    Perm v = w;
    std::swap(v.p[r - 1], v.p[s - 1]);
    // x_r * S_v
    {
        ExpMonoDict sub;
        trans_rec(v, n, sub, leaves, budget);
        for (auto& kv : sub) {
            ExpMono m = kv.first;
            if ((int)m.size() < r) m.resize(r, 0);
            m[r - 1] += 1;
            add_mono(out, std::move(m), kv.second);
        }
    }
    // + sum over i < r with v(i) < v(r), taking each value as a new maximum from the right
    int vr = v.p[r - 1], last = 0;
    for (int i = r - 1; i >= 1; --i) {
        int vi = v.p[i - 1];
        if (last < vi && vi < vr) {
            last = vi;
            Perm next = v;
            std::swap(next.p[i - 1], next.p[r - 1]);
            trans_rec(next, n, out, leaves, budget);
        }
    }
}

// Monk's rule: out += coeff * x_i * S_w summed over the dict, in S_n.
static void monk_add(int i, const std::unordered_map<Perm, int_coef_t, PermHash>& slc, int n, std::unordered_map<Perm, int_coef_t, PermHash>& out) {
    for (const auto& kv : slc) {
        const Perm& w = kv.first;
        int_coef_t c = kv.second;
        if (!c) continue;
        int wi = w.p[i - 1];
        int last = 0;
        for (int j = i - 1; j >= 1; --j) {
            int wj = w.p[j - 1];
            if (last < wj && wj < wi) {
                last = wj;
                Perm u = w;
                std::swap(u.p[j - 1], u.p[i - 1]);
                if ((out[u] -= c) == 0) out.erase(u);
            }
        }
        last = n + 1;
        for (int j = i + 1; j <= n; ++j) {
            int wj = w.p[j - 1];
            if (wi < wj && wj < last) {
                last = wj;
                Perm u = w;
                std::swap(u.p[i - 1], u.p[j - 1]);
                if ((out[u] += c) == 0) out.erase(u);
            }
        }
    }
}

// Horner over the Schubert basis: out += (sum of terms) * S_perm.
static void mult_poly_rec(std::vector<std::pair<ExpMono, int_coef_t>>& terms, int maxvar, const Perm& perm, int n, std::unordered_map<Perm, int_coef_t, PermHash>& out) {
    if (terms.empty()) return;
    if (maxvar == 0) {
        int_coef_t c = 0;
        for (auto& t : terms) c += t.second;
        if (c && (out[perm] += c) == 0) out.erase(perm);
        return;
    }
    std::vector<std::pair<ExpMono, int_coef_t>> lower, upper;
    int mv0 = 0, mv1 = 0;
    for (auto& t : terms) {
        if ((int)t.first.size() < maxvar) {
            mv0 = std::max(mv0, (int)t.first.size());
            lower.push_back(std::move(t));
        } else {
            t.first[maxvar - 1] -= 1;
            trim(t.first);
            mv1 = std::max(mv1, (int)t.first.size());
            upper.push_back(std::move(t));
        }
    }
    terms.clear();
    std::unordered_map<Perm, int_coef_t, PermHash> res1;
    mult_poly_rec(upper, mv1, perm, n, res1);
    monk_add(maxvar, res1, n, out);
    mult_poly_rec(lower, mv0, perm, n, out);
}

}  // namespace trans_detail

// Monomial expansion of S_w in S_n; throws TransBudgetExceeded once more than `budget` pipe dreams
// have been visited.
static ExpMonoDict trans_polynomial(const Perm& w, int n, long budget = LONG_MAX) {
    ExpMonoDict out;
    long leaves = 0;
    trans_detail::trans_rec(w, n, out, leaves, budget);
    return out;
}

namespace pd_detail {

typedef std::unordered_map<Perm, long, PermHash> CountMemo;

static long count_rec(const Perm& w, int n, CountMemo& memo, long budget);

// Row 1 of a pipe dream for w has crosses in columns 1..w(1)-1 and a subset of columns > w(1)
// whose decreasing word s_c (left multiplication, swapping values c, c+1) stays reduced against w.
// Enumerate the subsets for columns c, c-1, ..., m+1 and recurse on w with row 1 removed.
static long first_rows(Perm cur, int c, int m, int n, CountMemo& memo, long budget) {
    long total = 0;
    for (; c > m; --c) {
        int pc = 0, pc1 = 0;
        for (int i = 0; i < n; ++i) {
            if (cur.p[i] == c) pc = i;
            else if (cur.p[i] == c + 1) pc1 = i;
        }
        if (pc1 < pc) {
            Perm nxt = cur;
            std::swap(nxt.p[pc], nxt.p[pc1]);
            total += first_rows(nxt, c - 1, m, n, memo, budget);
            if (total > budget) throw TransBudgetExceeded{};
        }
    }
    // the forced crosses s_{m-1} ... s_1 move value m to 1; delete it and standardize
    Perm r = identity_perm();
    for (int i = 1; i < n; ++i) r.p[i - 1] = (uint8_t)(cur.p[i] > m ? cur.p[i] - 1 : cur.p[i]);
    return total + count_rec(r, n - 1, memo, budget);
}

static long count_rec(const Perm& w, int n, CountMemo& memo, long budget) {
    while (n > 1 && w.p[n - 1] == n) --n;
    if (n <= 1) return 1;
    auto it = memo.find(w);
    if (it != memo.end()) return it->second;
    long c = first_rows(w, n - 1, w.p[0], n, memo, budget);
    if (c > budget) throw TransBudgetExceeded{};
    memo.emplace(w, c);
    return c;
}

}  // namespace pd_detail

// S_w(1, ..., 1) = number of pipe dreams, by peeling off the first row (RCGraph.count_rc_graphs);
// memoized on the sub-permutations, so far cheaper than the transition expansion. Returns -1 as
// soon as any intermediate count exceeds `budget` (the root count dominates every sub-count).
static long pipe_dream_count(const Perm& w, int n, long budget = LONG_MAX) {
    pd_detail::CountMemo memo;
    try {
        return pd_detail::count_rec(w, n, memo, budget);
    } catch (const TransBudgetExceeded&) {
        return -1;
    }
}

// (sum of poly) * S_perm expanded in the Schubert basis of S_n.
static IntDict mult_poly_schubert(const ExpMonoDict& poly, const Perm& perm, int n) {
    std::vector<std::pair<ExpMono, int_coef_t>> terms(poly.begin(), poly.end());
    int maxvar = 0;
    for (auto& t : terms) maxvar = std::max(maxvar, (int)t.first.size());
    std::unordered_map<Perm, int_coef_t, PermHash> out;
    trans_detail::mult_poly_rec(terms, maxvar, perm, n, out);
    return IntDict(out.begin(), out.end());
}

// sum_u c_u S_u * S_v by expanding S_v into monomials (transition) and Horner-Monk into each S_u.
static IntDict schubmult_transition(const IntDict& perm_dict, const Perm& v, int n) {
    ExpMonoDict poly = trans_polynomial(v, n);
    std::unordered_map<Perm, int_coef_t, PermHash> acc;
    for (const auto& kv : perm_dict) {
        if (!kv.second) continue;
        for (const auto& t : mult_poly_schubert(poly, kv.first, n))
            if ((acc[t.first] += kv.second * t.second) == 0) acc.erase(t.first);
    }
    return IntDict(acc.begin(), acc.end());
}

// Number of v-path transitions in the setup for v: the size of the structure the v-path kernel
// iterates over (its cost scales like this to the ~0.8 power times the output size).
static long vpath_count(const MultSetup& S) {
    long c = 0;
    for (const auto& layer : S.vp.trans)
        for (const auto& steps : layer) c += (long)steps.size();
    return c;
}

// Cost-model dispatch for S_u * S_v.  Both kernels work on the factor e with fewer inversions
// (v-path: recurse on e; transition: expand e).  On exhaustive S_7 timings (log-ns, R^2 .90/.87)
//   transition ~ 7.70 + 0.74 log(1 + pd_e) + 0.50 log(1 + output)
//   v-path     ~ 7.91 + 0.83 log(1 + nvp_e) + 0.43 log(1 + output)
// with pd_e the pipe dreams of e and nvp_e the v-path transitions of e; on random, Grassmannian
// and many-descent pairs in S_8..S_12 the v-path cost grows faster, ~ nvp^1.0 against pd^0.6.
// The rule pd^0.6 < 0.25 nvp^1.1 is the compromise that minimizes the simulated total on both:
// S_7 1320 s (oracle 1179, v-path alone 1469, transition alone 2662), S_8..S_12 pairs 58 s
// (oracle 50, v-path 83, transition 426).  One v-path setup (reused by the kernel) and one
// pipe-dream count, which stops as soon as it passes the threshold, decide.
struct HybridChoice {
    bool transition;
    const Perm* expand;  // transition: factor expanded into monomials; v-path: factor recursed on
    const Perm* other;
    long pd, nvp;
};

// Largest pipe-dream count for which transition is predicted to beat a v-path run of size nvp
// (-1: never).
static long hybrid_pd_threshold(long nvp) {
    double lim = std::pow(0.247 * std::pow((double)(1 + nvp), 1.1), 1.0 / 0.6) - 1;
    return lim < 0 ? -1 : lim > 1e15 ? LONG_MAX : (long)lim;
}

// `setup` is left holding the v-path setup of the chosen factor.
static HybridChoice hybrid_choose(const Perm& u, const Perm& v, int n, MultSetup& setup) {
    bool on_u = perm_inv(u, n) < perm_inv(v, n);
    const Perm& e = on_u ? u : v;
    setup = mult_setup(e);
    long nvp = setup.trivial ? 0 : vpath_count(setup);
    long pd = pipe_dream_count(e, n, hybrid_pd_threshold(nvp));
    return HybridChoice{pd >= 0, &e, on_u ? &v : &u, pd, nvp};
}

// S_u * S_v by whichever of the two kernels the cost model predicts to be faster.
static IntDict schubmult_single_hybrid(const Perm& u, const Perm& v, int n) {
    MultSetup setup;
    HybridChoice c = hybrid_choose(u, v, n, setup);
    if (c.transition) return mult_poly_schubert(trans_polynomial(*c.expand, n), *c.other, n);
    return schubmult_single_with_setup({{*c.other, 1}}, setup, n);
}

// (sum_u c_u S_u) * S_v with cost-model dispatch.  A single term may work on either factor; for
// several terms only v can be expanded or recursed on.
static IntDict schubmult_hybrid(const IntDict& perm_dict, const Perm& v, int n) {
    const std::pair<Perm, int_coef_t>* only = nullptr;
    size_t nonzero = 0;
    for (const auto& kv : perm_dict)
        if (kv.second) { ++nonzero; only = &kv; }
    if (nonzero == 0) return {};
    if (nonzero == 1) {
        IntDict r = schubmult_single_hybrid(only->first, v, n);
        if (only->second != 1)
            for (auto& t : r) t.second *= only->second;
        return r;
    }
    MultSetup sv = mult_setup(v);
    long nvp = sv.trivial ? 0 : vpath_count(sv);
    if (pipe_dream_count(v, n, hybrid_pd_threshold(nvp)) >= 0) return schubmult_transition(perm_dict, v, n);
    return schubmult_single_with_setup(perm_dict, sv, n);
}
