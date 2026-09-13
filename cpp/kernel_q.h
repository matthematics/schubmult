// kernel_q.h: C++ port of schubmult_q_fast (library form; the CLI lives in schubmult_q_core.cpp).

#pragma once

#include "schub_quantum.h"

// ---------------------------------------------------------------------------
// q-polynomials with int coefficients
// ---------------------------------------------------------------------------

struct QPoly {
    std::vector<std::pair<Mono, int>> t;
    bool empty() const { return t.empty(); }
};

// dst += s * src * m
static void qpoly_addmul(QPoly& dst, const QPoly& src, int s, const Mono& m) {
    for (const auto& term : src.t) {
        Mono nm = mono_mul(term.first, m);
        int c = s * term.second;
        bool found = false;
        for (auto& d : dst.t)
            if (d.first == nm) {
                d.second += c;
                found = true;
                break;
            }
        if (!found) dst.t.push_back({nm, c});
    }
}

static void qpoly_prune(QPoly& p) {
    p.t.erase(std::remove_if(p.t.begin(), p.t.end(), [](const std::pair<Mono, int>& x) { return x.second == 0; }), p.t.end());
}

static std::string qpoly_str(QPoly p) {
    std::sort(p.t.begin(), p.t.end());
    std::string s;
    for (size_t idx = 0; idx < p.t.size(); ++idx) {
        int c = p.t[idx].second;
        const Mono& m = p.t[idx].first;
        std::string mono;
        for (int i = 0; i < MAXN; ++i) {
            if (!m.e[i]) continue;
            if (!mono.empty()) mono += "*";
            mono += "q_" + std::to_string(i + 1);
            if (m.e[i] > 1) mono += "**" + std::to_string((int)m.e[i]);
        }
        int ac = c < 0 ? -c : c;
        std::string term;
        if (mono.empty())
            term = std::to_string(ac);
        else if (ac == 1)
            term = mono;
        else
            term = std::to_string(ac) + "*" + mono;
        if (idx == 0)
            s += (c < 0 ? "-" : "") + term;
        else
            s += (c < 0 ? " - " : " + ") + term;
    }
    return s;
}

// ---------------------------------------------------------------------------
// schubmult_q_fast
// ---------------------------------------------------------------------------

typedef std::vector<std::pair<Perm, QPoly>> QDict;
typedef PermTable<QPoly> QTable;
typedef QTable::VSum QVSum;

static QPoly& local_get_or_insert(std::vector<QVSum>& vec, uint32_t v) {
    for (QVSum& e : vec)
        if (e.first == v) return e.second;
    vec.push_back({v, QPoly()});
    return vec.back().second;
}

static QDict schubmult_q_fast(const QDict& perm_dict, const Perm& v) {
    if (perm_inv(v, MAXN) == 0) return perm_dict;
    MultSetup S = mult_setup(v, /*medium=*/true);
    if (S.trivial) return perm_dict;
    const VPaths& vp = S.vp;
    const std::vector<int>& th = S.th;
    int thL = (int)th.size();
    if (vp.id_start == UINT32_MAX) return {};
    uint32_t id_vmu = 0;  // level[thL] == {vmu}
    const Mono one = mono_one();

    std::unordered_map<Perm, QPoly, PermHash> result;
    QTable tabA, tabB;
    QTable* A = &tabA;
    QTable* B = &tabB;

    std::vector<QUp> newperms, keys;
    std::vector<std::vector<QUp>> second;
    std::vector<std::vector<QVSum>> nps0;  // newpathsums0: per key, sums over level[index+1]

    for (const auto& kv : perm_dict) {
        const Perm& u = kv.first;
        if (kv.second.empty()) continue;
        int inv_u = perm_inv(u, MAXN);

        A->reset();
        A->sums[A->intern(u, inv_u)].push_back({vp.id_start, kv.second});

        for (int index = 0; index < thL; ++index) {
            if (index > 0 && th[index - 1] == th[index]) continue;  // consumed by the double step below
            int k = th[index];
#ifdef STATS
            double t0 = now_s();
#endif
            B->reset();
            if (index < thL - 1 && th[index] == th[index + 1]) {
                int mx0 = vp.mx_th[index], mx1 = vp.mx_th[index + 1];
                const auto& trans0 = vp.trans[index];
                const auto& trans1 = vp.trans[index + 1];
                for (uint32_t sid = 0; sid < A->count; ++sid) {
                    const std::vector<QVSum>& sums = A->sums[sid];
                    if (sums.empty()) continue;
                    const Perm& up = A->perms[sid];
                    double_elem_sym_q(up, mx0, mx1, k, keys, second);
                    nps0.assign(keys.size(), {});
                    for (const QVSum& sv : sums)
                        for (const Trans& tr : trans0[sv.first])
                            for (size_t kk = 0; kk < keys.size(); ++kk)
                                if (keys[kk].udiff + tr.vdiff == k) qpoly_addmul(local_get_or_insert(nps0[kk], tr.v2), sv.second, tr.s, keys[kk].q);
                    for (size_t kk = 0; kk < keys.size(); ++kk) {
                        if (nps0[kk].empty()) continue;
                        for (const QUp& e2 : second[kk]) {
                            long id = -1;
                            for (const QVSum& sv : nps0[kk]) {
                                if (sv.second.empty()) continue;
                                for (const Trans& tr : trans1[sv.first]) {
                                    if (e2.udiff + tr.vdiff != th[index + 1]) continue;
                                    if (id < 0) id = B->intern(e2.perm, perm_inv(e2.perm, MAXN));
                                    qpoly_addmul(B->get_or_insert((uint32_t)id, tr.v2), sv.second, tr.s, e2.q);
                                }
                            }
                        }
                    }
                }
            } else {
                const auto& trans = vp.trans[index];
                for (uint32_t sid = 0; sid < A->count; ++sid) {
                    const std::vector<QVSum>& sums = A->sums[sid];
                    if (sums.empty()) continue;
                    const Perm& up = A->perms[sid];
                    int inv_up = A->inv[sid];
                    int p = std::min(vp.mx_th[index], (S.inv_mu - (inv_up - inv_u)) - S.inv_vmu);
                    elem_sym_perms_q(up, p, k, newperms);
                    for (const QUp& e : newperms) {
                        long id = -1;
                        for (const QVSum& sv : sums) {
                            for (const Trans& tr : trans[sv.first]) {
                                if (e.udiff + tr.vdiff != k) continue;
                                if (id < 0) id = B->intern(e.perm, perm_inv(e.perm, MAXN));
                                qpoly_addmul(B->get_or_insert((uint32_t)id, tr.v2), sv.second, tr.s, e.q);
                            }
                        }
                    }
                }
            }

            // drop cancelled terms and entries
            size_t alive = 0;
            for (auto& vec : B->sums) {
                for (QVSum& e : vec) qpoly_prune(e.second);
                vec.erase(std::remove_if(vec.begin(), vec.end(), [](const QVSum& e) { return e.second.empty(); }), vec.end());
                if (!vec.empty()) ++alive;
            }
            for (uint32_t sid = 0; sid < B->count; ++sid)
                if (B->perms[sid].p[MAXN - 1] != MAXN) die("permutation reached length MAXN; rebuild with a larger -DMAXN");
            std::swap(A, B);
#ifdef STATS
            std::fprintf(stderr, "idx %2d k=%2d targets=%8u states_out=%8zu  %.3fs\n", index, k, A->count, alive, now_s() - t0);
#endif
        }

        for (uint32_t sid = 0; sid < A->count; ++sid)
            if (QPoly* c = A->find(sid, id_vmu)) qpoly_addmul(result[A->perms[sid]], *c, 1, one);
    }

    QDict out;
    out.reserve(result.size());
    for (auto& kv : result) {
        qpoly_prune(kv.second);
        if (!kv.second.empty()) out.push_back(kv);
    }
    std::sort(out.begin(), out.end(), [](const std::pair<Perm, QPoly>& x, const std::pair<Perm, QPoly>& y) { return x.first < y.first; });
    return out;
}

