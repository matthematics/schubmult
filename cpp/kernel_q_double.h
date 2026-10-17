// kernel_q_double.h: C++ port of schubmult_q_double_fast (library form; the CLI lives in schubmult_q_double_core.cpp).

#pragma once

#include "schub_quantum.h"
#include "schub_symbolic.h"

#ifndef SCHUB_EXPRDICT
#define SCHUB_EXPRDICT
typedef std::vector<std::pair<Perm, Expr>> ExprDict;
#endif

static Expr mono_expr(const Mono& m, const ElemSymCache& esc) {
    ExprVec fs;
    for (int i = 0; i < MAXN; ++i)
        if (m.e[i]) fs.push_back(m.e[i] == 1 ? ElemSymCache::need(esc.Q, i + 1) : Expr(ex_pow(ElemSymCache::need(esc.Q, i + 1), m.e[i])));
    return ex_mul(fs);
}

// Shadow of a q-monomial from its exponents: mono_expr builds a fresh object each time, so it cannot
// be looked up by address; the q_i themselves are stable (esc.Q) and are.
static ShadowVal mono_shadow(const Mono& m, const ElemSymCache& esc, Shadow& shadow) {
    ShadowVal r;
    for (int i = 0; i < SHADOW_POINTS; ++i) r.v[i] = 1;
    for (int i = 0; i < MAXN; ++i)
        if (m.e[i]) {
            const ShadowVal& q = shadow.of(ElemSymCache::need(esc.Q, i + 1));
            for (int k = 0; k < m.e[i]; ++k) r = r * q;
        }
    return r;
}

// elem_sym_func_q(k, i, u1, u2, v1, v2, udiff, vdiff): nullptr when it is 0, ex_one() when 1.
static Expr elem_sym_func_q(int k, const Perm& u1, const Perm& u2, int udiff, const Trans& tr, const std::vector<int>& zs, ElemSymCache& esc, std::vector<int>& yidx) {
    int newk = k - udiff;
    if (newk < tr.vdiff) return Expr();
    if (newk == tr.vdiff) return ex_one();
    yvars_of(u1, u2, k, yidx);
    return esc.get(newk - tr.vdiff, newk, yidx, zs);
}

// Shadow of an elem_sym_func_q value: the cached e_p's are stable objects; ex_one() is handled by hand.
static ShadowVal esf_shadow(const Expr& esf, int k, int udiff, const Trans& tr, Shadow& shadow) {
    if (k - udiff == tr.vdiff) {
        ShadowVal one;
        for (int i = 0; i < SHADOW_POINTS; ++i) one.v[i] = 1;
        return one;
    }
    return shadow.of(esf);
}

// ---------------------------------------------------------------------------
// A partial sum, split by the q-monomials introduced by this multiplication: for each monomial the
// pending terms of its (y, z)-coefficient and that coefficient's shadow. Input coefficients may
// themselves contain q (from earlier multiplications); those q's stay inside the coefficients, under
// the monomial 1. Keeping the q-monomials apart lets a component that has become zero be dropped
// while the others go on, and the output is assembled as sum_m q^m c_m at the end.
// ---------------------------------------------------------------------------

struct QPending {
    ExprVec terms;
    ShadowVal val;
};

struct QExpr {
    std::vector<std::pair<Mono, QPending>> t;
    bool empty() const { return t.empty(); }
};

static void qexpr_push(QExpr& dst, const Mono& m, const Expr& term, const ShadowVal& tv) {
    for (auto& d : dst.t)
        if (d.first == m) {
            d.second.terms.push_back(term);
            d.second.val += tv;
            return;
        }
    dst.t.push_back({m, QPending{ExprVec{term}, tv}});
}

// Sum each component's pending terms; drop the structurally zero ones and, with a shadow, those whose
// coefficient vanishes at the sample points. If the whole sum vanishes there (possible only through
// q's carried inside the coefficients), drop everything.
static void qexpr_collapse(QExpr& qe, const ElemSymCache& esc, Shadow* shadow) {
    ShadowVal total;
    for (auto& d : qe.t) {
        QPending& pd = d.second;
        Expr s = ex_add(pd.terms);
        pd.terms.clear();
        if (is_zero(s) || (shadow && pd.val.is_zero())) continue;
        pd.terms.push_back(s);
        if (shadow) total += mono_shadow(d.first, esc, *shadow) * pd.val;
    }
    qe.t.erase(std::remove_if(qe.t.begin(), qe.t.end(), [](const std::pair<Mono, QPending>& d) { return d.second.terms.empty(); }), qe.t.end());
    if (shadow && !qe.t.empty() && total.is_zero()) qe.t.clear();
}

static Expr qexpr_value(const QExpr& qe, const ElemSymCache& esc) {
    ExprVec parts;
    for (const auto& d : qe.t) {
        const Expr& c = d.second.terms[0];
        bool one = true;
        for (int i = 0; i < MAXN; ++i)
            if (d.first.e[i]) one = false;
        parts.push_back(one ? c : ex_mul(c, mono_expr(d.first, esc)));
    }
    return ex_add(parts);
}

static ExprDict schubmult_q_double_fast(const ExprDict& perm_dict, const Perm& v, ElemSymCache& esc, Shadow* shadow = nullptr) {
    if (perm_inv(v, MAXN) == 0) return perm_dict;
    MultSetup S = mult_setup(v, /*medium=*/true);
    if (S.trivial) return perm_dict;
    const VPaths& vp = S.vp;
    const std::vector<int>& th = S.th;
    int thL = (int)th.size();
    if (vp.id_start == UINT32_MAX) return {};
    uint32_t id_vmu = 0;  // level[thL] == {vmu}
    ZIdx zidx = compute_zidx(vp);

    typedef PermTable<QExpr> Table;
    typedef Table::VSum VSum;
    std::unordered_map<Perm, QExpr, PermHash> result;

    Table tabA, tabB;
    Table* A = &tabA;
    Table* B = &tabB;

    std::vector<QUp> newperms, keys;
    std::vector<std::vector<QUp>> second;
    std::vector<std::vector<VSum>> nps0;  // newpathsums0: per key, pending sums over level[index+1]
    std::vector<int> yidx;

    auto local_get = [](std::vector<VSum>& vec, uint32_t v2) -> QExpr& {
        for (VSum& e : vec)
            if (e.first == v2) return e.second;
        vec.push_back({v2, QExpr()});
        return vec.back().second;
    };
    auto collapse_all = [&](std::vector<VSum>& vec) {
        for (VSum& e : vec) qexpr_collapse(e.second, esc, shadow);
        vec.erase(std::remove_if(vec.begin(), vec.end(), [](const VSum& e) { return e.second.empty(); }), vec.end());
    };
    // one term for each q-component of a source state: coefficient * esf, with sign, into dst under the
    // component's monomial times qm
    auto spread = [&](const QExpr& src, const Expr& esf, const ShadowVal& esfv, bool esf_one, int sign, const Mono& qm, QExpr& dst) {
        for (const auto& d : src.t) {
            const Expr& sumval = d.second.terms[0];
            Expr term = esf_one ? sumval : ex_mul(sumval, esf);
            ShadowVal tv;
            if (shadow) tv = esf_one ? d.second.val : d.second.val * esfv;
            if (sign < 0) {
                term = ex_neg(term);
                tv = -tv;
            }
            qexpr_push(dst, mono_mul(d.first, qm), term, tv);
        }
    };

    for (const auto& kv : perm_dict) {
        const Perm& u = kv.first;
        if (is_zero(kv.second)) continue;
        int inv_u = perm_inv(u, MAXN);

        QExpr start;
        start.t.push_back({mono_one(), QPending{ExprVec{kv.second}, ShadowVal()}});
        if (shadow) {
            start.t[0].second.val = shadow->of(kv.second);
            if (start.t[0].second.val.is_zero()) continue;
        }
        A->reset();
        A->sums[A->intern(u, inv_u)].push_back({vp.id_start, std::move(start)});

        for (int index = 0; index < thL; ++index) {
            if (index > 0 && th[index - 1] == th[index]) continue;  // consumed by the double step below
            int k = th[index];
#ifdef STATS
            double t0 = now_s();
#endif
            B->reset();
            if (index < thL - 1 && th[index] == th[index + 1]) {
                int k1 = th[index + 1];
                const auto& trans0 = vp.trans[index];
                const auto& trans1 = vp.trans[index + 1];
                for (uint32_t sid = 0; sid < A->count; ++sid) {
                    const std::vector<VSum>& sums = A->sums[sid];
                    if (sums.empty()) continue;
                    const Perm& up = A->perms[sid];
                    double_elem_sym_q(up, vp.mx_th[index], vp.mx_th[index + 1], k, keys, second);
                    nps0.assign(keys.size(), {});
                    for (const VSum& sv : sums) {
                        const auto& trs = trans0[sv.first];
                        for (size_t t = 0; t < trs.size(); ++t) {
                            const Trans& tr = trs[t];
                            for (size_t kk = 0; kk < keys.size(); ++kk) {
                                Expr esf = elem_sym_func_q(k, up, keys[kk].perm, keys[kk].udiff, tr, zidx[index][sv.first][t], esc, yidx);
                                if (esf.is_null()) continue;
                                bool one = k - keys[kk].udiff == tr.vdiff;
                                ShadowVal esfv;
                                if (shadow) esfv = esf_shadow(esf, k, keys[kk].udiff, tr, *shadow);
                                spread(sv.second, esf, esfv, one, tr.s, keys[kk].q, local_get(nps0[kk], tr.v2));
                            }
                        }
                    }
                    for (size_t kk = 0; kk < keys.size(); ++kk) {
                        collapse_all(nps0[kk]);
                        if (nps0[kk].empty()) continue;
                        const Perm& up1 = keys[kk].perm;
                        for (const VSum& sv : nps0[kk]) {
                            const auto& trs = trans1[sv.first];
                            for (size_t t = 0; t < trs.size(); ++t) {
                                const Trans& tr = trs[t];
                                for (const QUp& e2 : second[kk]) {
                                    Expr esf = elem_sym_func_q(k1, up1, e2.perm, e2.udiff, tr, zidx[index + 1][sv.first][t], esc, yidx);
                                    if (esf.is_null()) continue;
                                    bool one = k1 - e2.udiff == tr.vdiff;
                                    ShadowVal esfv;
                                    if (shadow) esfv = esf_shadow(esf, k1, e2.udiff, tr, *shadow);
                                    uint32_t id = B->intern(e2.perm, perm_inv(e2.perm, MAXN));
                                    spread(sv.second, esf, esfv, one, tr.s, e2.q, B->get_or_insert(id, tr.v2));
                                }
                            }
                        }
                    }
                }
            } else {
                const auto& trans = vp.trans[index];
                for (uint32_t sid = 0; sid < A->count; ++sid) {
                    const std::vector<VSum>& sums = A->sums[sid];
                    if (sums.empty()) continue;
                    const Perm& up = A->perms[sid];
                    int inv_up = A->inv[sid];
                    int p = std::min(vp.mx_th[index], (S.inv_mu - (inv_up - inv_u)) - S.inv_vmu);
                    elem_sym_perms_q(up, p, k, newperms);
                    for (const QUp& e : newperms) {
                        long id = -1;
                        for (const VSum& sv : sums) {
                            const auto& trs = trans[sv.first];
                            for (size_t t = 0; t < trs.size(); ++t) {
                                const Trans& tr = trs[t];
                                Expr esf = elem_sym_func_q(k, up, e.perm, e.udiff, tr, zidx[index][sv.first][t], esc, yidx);
                                if (esf.is_null()) continue;
                                bool one = k - e.udiff == tr.vdiff;
                                ShadowVal esfv;
                                if (shadow) esfv = esf_shadow(esf, k, e.udiff, tr, *shadow);
                                if (id < 0) id = B->intern(e.perm, perm_inv(e.perm, MAXN));
                                spread(sv.second, esf, esfv, one, tr.s, e.q, B->get_or_insert((uint32_t)id, tr.v2));
                            }
                        }
                    }
                }
            }

            size_t alive = 0;
            for (auto& vec : B->sums) {
                collapse_all(vec);
                if (!vec.empty()) ++alive;
            }
            for (uint32_t sid = 0; sid < B->count; ++sid)
                if (B->perms[sid].p[MAXN - 1] != MAXN) die("permutation reached length MAXN; rebuild with a larger -DMAXN");
            std::swap(A, B);
#ifdef STATS
            std::fprintf(stderr, "idx %2d k=%2d targets=%8u states_out=%8zu esf_cache=%zu  %.3fs\n", index, k, A->count, alive, esc.memo.size(), now_s() - t0);
#endif
        }

        for (uint32_t sid = 0; sid < A->count; ++sid)
            if (QExpr* c = A->find(sid, id_vmu)) {
                QExpr& r = result[A->perms[sid]];
                for (const auto& d : c->t) qexpr_push(r, d.first, d.second.terms[0], d.second.val);
            }
    }

    ExprDict out;
    out.reserve(result.size());
    for (auto& kv : result) {
        qexpr_collapse(kv.second, esc, shadow);
        if (kv.second.empty()) continue;
        out.push_back({kv.first, qexpr_value(kv.second, esc)});
    }
    std::sort(out.begin(), out.end(), [](const std::pair<Perm, Expr>& x, const std::pair<Perm, Expr>& y) { return x.first < y.first; });
    return out;
}
