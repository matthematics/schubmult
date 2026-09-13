// kernel_single.h: C++ port of schubmult_py (library form; the CLI lives in schubmult_core.cpp).

#pragma once

#include "schub_common.h"

typedef int int_coef_t;
typedef std::vector<std::pair<Perm, int_coef_t>> IntDict;

// n: every permutation in the answer lies in S_n (the caller derives it from input sizes);
// paths leaving S_n can never come back, so positions >= n are never touched.
static IntDict schubmult_single(const IntDict& perm_dict, const Perm& v, int n) {
    MultSetup S = mult_setup(v);
    if (S.trivial) return perm_dict;
    const VPaths& vp = S.vp;
    const std::vector<int>& th = S.th;
    int thL = (int)th.size();
    if (vp.id_start == UINT32_MAX) return {};
    uint32_t id_vmu = 0;  // level[thL] == {vmu}

    typedef PermTable<int_coef_t> Table;
    typedef Table::VSum VSum;
    std::unordered_map<Perm, int_coef_t, PermHash> result;

    // Two tables alternate: A holds the current step's permutations, B receives the next.
    Table tabA, tabB;
    Table* A = &tabA;
    Table* B = &tabB;

    // scratch
    std::vector<std::pair<Perm, int>> newperms;
    std::vector<int_coef_t> loc;
    std::vector<uint32_t> loc_touched;
    std::vector<std::vector<VSum>> contrib;  // contrib[d] = nonzero (v2, c) pairs of this state at udiff d

    for (const auto& kv : perm_dict) {
        const Perm& u = kv.first;
        int_coef_t val = kv.second;
        if (val == 0) continue;
        int inv_u = perm_inv(u, n);

        A->reset();
        A->sums[A->intern(u, inv_u)].push_back({vp.id_start, val});

        for (int index = 0; index < thL; ++index) {
            int k = th[index];
            const auto& trans = vp.trans[index];
            uint32_t nlev = (uint32_t)vp.level[index + 1].size();
#ifdef STATS
            double t0 = now_s(), t_gen = 0;
            size_t n_newperms = 0, n_adds = 0;
#endif
            B->reset();
            if ((int)contrib.size() < k + 1) contrib.resize(k + 1);
            loc.assign((size_t)(vp.mx_th[index] + 1) * nlev, 0);
            loc_touched.clear();

            for (uint32_t sid = 0; sid < A->count; ++sid) {
                const std::vector<VSum>& sums = A->sums[sid];
                if (sums.empty()) continue;
                const Perm& up = A->perms[sid];
                int inv_up = A->inv[sid];
                int p = std::min(vp.mx_th[index], S.inv_mu - S.inv_vmu - (inv_up - inv_u));
                if (p < 0) p = 0;

                // aggregate this state's contributions by (udiff, v2) first; udiffs with none
                // need no permutations generated at all
                for (const VSum& sv : sums) {
                    for (const Trans& tr : trans[sv.first]) {
                        int d = k - tr.vdiff;
                        if (d < 0 || d > p) continue;
                        uint32_t idx = (uint32_t)d * nlev + tr.v2;
                        if (loc[idx] == 0) loc_touched.push_back(idx);
                        loc[idx] += tr.s * sv.second;
                    }
                }
                for (int d = 0; d <= p; ++d) contrib[d].clear();
                int pmax = -1;
                for (uint32_t idx : loc_touched) {
                    int_coef_t c = loc[idx];
                    loc[idx] = 0;
                    if (c == 0) continue;
                    uint32_t d = idx / nlev, v2 = idx - d * nlev;
                    contrib[d].push_back({v2, c});
                    pmax = std::max(pmax, (int)d);
                }
                loc_touched.clear();
                if (pmax < 0) continue;
#ifdef STATS
                double tg = now_s();
#endif
                elem_sym_perms(up, pmax, k, n, newperms);
#ifdef STATS
                t_gen += now_s() - tg;
                n_newperms += newperms.size();
#endif
                for (const auto& npd : newperms) {
                    const std::vector<VSum>& cs = contrib[npd.second];
                    if (cs.empty()) continue;
                    uint32_t id = B->intern(npd.first, inv_up + npd.second);
                    for (const VSum& c : cs) B->get_or_insert(id, c.first) += c.second;
#ifdef STATS
                    n_adds += cs.size();
#endif
                }
            }

            // drop cancelled entries (permutations with no surviving sums are skipped next step)
            size_t alive = 0, live_sums = 0;
            for (auto& vec : B->sums) {
                vec.erase(std::remove_if(vec.begin(), vec.end(), [](const VSum& e) { return e.second == 0; }), vec.end());
                if (!vec.empty()) {
                    ++alive;
                    live_sums += vec.size();
                }
            }
            std::swap(A, B);
#ifdef STATS
            std::fprintf(stderr, "idx %2d k=%2d |L|=%5u newperms=%10zu targets=%8u adds=%11zu states_out=%8zu sums_out=%10zu  %.3fs (gen %.3fs)\n", index, k, nlev, n_newperms, A->count, n_adds,
                         alive, live_sums, now_s() - t0, t_gen);
#endif
        }

        for (uint32_t sid = 0; sid < A->count; ++sid)
            if (int_coef_t* c = A->find(sid, id_vmu)) result[A->perms[sid]] += *c;
    }

    IntDict out;
    out.reserve(result.size());
    for (auto& kv : result)
        if (kv.second != 0) out.push_back(kv);
    std::sort(out.begin(), out.end(), [](const std::pair<Perm, int_coef_t>& x, const std::pair<Perm, int_coef_t>& y) { return x.first < y.first; });
    return out;
}

