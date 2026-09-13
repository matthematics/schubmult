// kernel_double.h: C++ port of schubmult_double (library form; the CLI lives in schubmult_double_core.cpp).

#pragma once

#include "schub_symbolic.h"

// ---------------------------------------------------------------------------
// schubmult_double
// ---------------------------------------------------------------------------

#ifndef SCHUB_EXPRDICT
#define SCHUB_EXPRDICT
typedef std::vector<std::pair<Perm, Expr>> ExprDict;
#endif

static ExprDict schubmult_double(const ExprDict& perm_dict, const Perm& v, int n, ElemSymCache& esc) {
    MultSetup S = mult_setup(v);
    if (S.trivial) return perm_dict;
    const VPaths& vp = S.vp;
    const std::vector<int>& th = S.th;
    int thL = (int)th.size();
    if (vp.id_start == UINT32_MAX) return {};
    uint32_t id_vmu = 0;  // level[thL] == {vmu}

    ZIdx zidx = compute_zidx(vp);

    // Each sum entry holds the pending terms; at the end of a level they are summed
    // (SymEngine's Add merges like terms at the top level, as Python's + does) and
    // structurally zero entries are dropped. Nothing is ever expanded.
    typedef PermTable<ExprVec> Table;
    typedef Table::VSum VSum;
    std::unordered_map<Perm, ExprVec, PermHash> result;

    Table tabA, tabB;
    Table* A = &tabA;
    Table* B = &tabB;

    std::vector<std::pair<Perm, int>> newperms;
    std::vector<int> yidx;

    for (const auto& kv : perm_dict) {
        const Perm& u = kv.first;
        if (is_zero(kv.second)) continue;
        int inv_u = perm_inv(u, n);

        A->reset();
        A->sums[A->intern(u, inv_u)].push_back({vp.id_start, ExprVec{kv.second}});

        for (int index = 0; index < thL; ++index) {
            int k = th[index];
            const auto& trans = vp.trans[index];
#ifdef STATS
            double t0 = now_s();
            size_t n_terms = 0;
#endif
            B->reset();

            for (uint32_t sid = 0; sid < A->count; ++sid) {
                const std::vector<VSum>& sums = A->sums[sid];
                if (sums.empty()) continue;
                const Perm& up = A->perms[sid];
                int inv_up = A->inv[sid];
                int p = std::min(vp.mx_th[index], S.inv_mu - S.inv_vmu - (inv_up - inv_u));
                if (p < 0) p = 0;
                elem_sym_perms(up, p, k, n, newperms);

                for (const auto& npd : newperms) {
                    const Perm& up2 = npd.first;
                    int udiff = npd.second, newk = k - udiff;
                    yvars_of(up, up2, k, yidx);
                    long id = -1;
                    for (const VSum& sv : sums) {
                        const Expr& sumval = sv.second[0];
                        const auto& trs = trans[sv.first];
                        for (size_t t = 0; t < trs.size(); ++t) {
                            const Trans& tr = trs[t];
                            if (newk < tr.vdiff) continue;  // elem_sym_func == 0
                            Expr term;
                            if (newk == tr.vdiff) {
                                term = sumval;  // elem_sym_func == 1
                            } else {
                                const Expr& esf = esc.get(newk - tr.vdiff, newk, yidx, zidx[index][sv.first][t]);
                                if (is_zero(esf)) continue;
                                term = ex_mul(sumval, esf);
                            }
                            if (tr.s < 0) term = ex_neg(term);
                            if (id < 0) id = B->intern(up2, inv_up + udiff);
                            B->get_or_insert((uint32_t)id, tr.v2).push_back(term);
#ifdef STATS
                            ++n_terms;
#endif
                        }
                    }
                }
            }

            // sum each entry's pending terms; drop the zeros
            size_t alive = 0, live_sums = 0;
            for (auto& vec : B->sums) {
                for (VSum& e : vec) {
                    Expr s = ex_add(e.second);
                    e.second.clear();
                    if (!is_zero(s)) e.second.push_back(s);
                }
                vec.erase(std::remove_if(vec.begin(), vec.end(), [](const VSum& e) { return e.second.empty(); }), vec.end());
                if (!vec.empty()) {
                    ++alive;
                    live_sums += vec.size();
                }
            }
            std::swap(A, B);
#ifdef STATS
            std::fprintf(stderr, "idx %2d k=%2d terms=%10zu targets=%8u states_out=%8zu sums_out=%10zu esf_cache=%zu  %.3fs\n", index, k, n_terms, A->count, alive, live_sums, esc.memo.size(),
                         now_s() - t0);
#endif
        }

        for (uint32_t sid = 0; sid < A->count; ++sid)
            if (ExprVec* c = A->find(sid, id_vmu)) result[A->perms[sid]].push_back((*c)[0]);
    }

    ExprDict out;
    out.reserve(result.size());
    for (auto& kv : result) {
        Expr s = ex_add(kv.second);
        if (!is_zero(s)) out.push_back({kv.first, s});
    }
    std::sort(out.begin(), out.end(), [](const std::pair<Perm, Expr>& x, const std::pair<Perm, Expr>& y) { return x.first < y.first; });
    return out;
}

