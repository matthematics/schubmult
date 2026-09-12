// schubmult_q_double_core: C++ port of schubmult.mult.quantum_double.schubmult_q_double_fast
// (products of quantum double Schubert polynomials; no parabolic / nil-Hecke / positivity yet).
//
// Build:  make build/schubmult_q_double_core     (needs SymEngine C++ headers/lib, e.g. from conda)
// Usage:  ./build/schubmult_q_double_core 3 1 2 - 2 1 3                (z = y, the default)
//         ./build/schubmult_q_double_core --mixed-var 3 1 2 - 2 1 3    (S_u(x; y) * S_v(x; z))
//         ./build/schubmult_q_double_core --code 2 0 - 1 0
// Output lines are "(w1, w2, ...)  polynomial", sorted by (inv(w), w), zero coefficients dropped.

#include "schub_quantum.h"
#include "schub_symbolic.h"

typedef std::vector<std::pair<Perm, Expr>> PermDict;

static Expr mono_expr(const Mono& m, const ElemSymCache& esc) {
    vec_basic fs;
    for (int i = 0; i < MAXN; ++i)
        if (m.e[i]) fs.push_back(m.e[i] == 1 ? esc.Q[i + 1] : Expr(SymEngine::pow(esc.Q[i + 1], SymEngine::integer(m.e[i]))));
    return SymEngine::mul(fs);
}

// elem_sym_func_q(k, i, u1, u2, v1, v2, udiff, vdiff): nullptr when it is 0, SymEngine::one when 1.
static Expr elem_sym_func_q(int k, const Perm& u1, const Perm& u2, int udiff, const Trans& tr, const std::vector<int>& zs, ElemSymCache& esc, std::vector<int>& yidx) {
    int newk = k - udiff;
    if (newk < tr.vdiff) return Expr();
    if (newk == tr.vdiff) return SymEngine::one;
    yvars_of(u1, u2, k, yidx);
    return esc.get(newk - tr.vdiff, newk, yidx, zs);
}

static PermDict schubmult_q_double_fast(const PermDict& perm_dict, const Perm& v, ElemSymCache& esc) {
    if (perm_inv(v, MAXN) == 0) return perm_dict;
    MultSetup S = mult_setup(v, /*medium=*/true);
    if (S.trivial) return perm_dict;
    const VPaths& vp = S.vp;
    const std::vector<int>& th = S.th;
    int thL = (int)th.size();
    if (vp.id_start == UINT32_MAX) return {};
    uint32_t id_vmu = 0;  // level[thL] == {vmu}
    ZIdx zidx = compute_zidx(vp);

    // Each sum entry holds the pending terms; at the end of a level they are summed and
    // structurally zero entries are dropped. Nothing is ever expanded.
    typedef PermTable<vec_basic> Table;
    typedef Table::VSum VSum;
    std::unordered_map<Perm, vec_basic, PermHash> result;

    Table tabA, tabB;
    Table* A = &tabA;
    Table* B = &tabB;

    std::vector<QUp> newperms, keys;
    std::vector<std::vector<QUp>> second;
    std::vector<Expr> keyq, secq;
    std::vector<std::vector<VSum>> nps0;  // newpathsums0: per key, pending sums over level[index+1]
    std::vector<int> yidx;

    auto local_push = [](std::vector<VSum>& vec, uint32_t v2, const Expr& term) {
        for (VSum& e : vec)
            if (e.first == v2) {
                e.second.push_back(term);
                return;
            }
        vec.push_back({v2, vec_basic{term}});
    };
    auto collapse = [](std::vector<VSum>& vec) {
        for (VSum& e : vec) {
            Expr s = SymEngine::add(e.second);
            e.second.clear();
            if (!is_zero(s)) e.second.push_back(s);
        }
        vec.erase(std::remove_if(vec.begin(), vec.end(), [](const VSum& e) { return e.second.empty(); }), vec.end());
    };

    for (const auto& kv : perm_dict) {
        const Perm& u = kv.first;
        if (is_zero(kv.second)) continue;
        int inv_u = perm_inv(u, MAXN);

        A->reset();
        A->sums[A->intern(u, inv_u)].push_back({vp.id_start, vec_basic{kv.second}});

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
                    keyq.clear();
                    for (const QUp& e1 : keys) keyq.push_back(mono_expr(e1.q, esc));
                    nps0.assign(keys.size(), {});
                    for (const VSum& sv : sums) {
                        const Expr& sumval = sv.second[0];
                        const auto& trs = trans0[sv.first];
                        for (size_t t = 0; t < trs.size(); ++t) {
                            const Trans& tr = trs[t];
                            for (size_t kk = 0; kk < keys.size(); ++kk) {
                                Expr esf = elem_sym_func_q(k, up, keys[kk].perm, keys[kk].udiff, tr, zidx[index][sv.first][t], esc, yidx);
                                if (esf.is_null()) continue;
                                Expr term = SymEngine::mul(vec_basic{sumval, esf, keyq[kk]});
                                if (tr.s < 0) term = SymEngine::neg(term);
                                local_push(nps0[kk], tr.v2, term);
                            }
                        }
                    }
                    for (size_t kk = 0; kk < keys.size(); ++kk) {
                        collapse(nps0[kk]);
                        if (nps0[kk].empty()) continue;
                        const Perm& up1 = keys[kk].perm;
                        secq.clear();
                        for (const QUp& e2 : second[kk]) secq.push_back(mono_expr(e2.q, esc));
                        for (const VSum& sv : nps0[kk]) {
                            const Expr& sumval = sv.second[0];
                            const auto& trs = trans1[sv.first];
                            for (size_t t = 0; t < trs.size(); ++t) {
                                const Trans& tr = trs[t];
                                for (size_t s2 = 0; s2 < second[kk].size(); ++s2) {
                                    const QUp& e2 = second[kk][s2];
                                    Expr esf = elem_sym_func_q(k1, up1, e2.perm, e2.udiff, tr, zidx[index + 1][sv.first][t], esc, yidx);
                                    if (esf.is_null()) continue;
                                    Expr term = SymEngine::mul(vec_basic{sumval, esf, secq[s2]});
                                    if (tr.s < 0) term = SymEngine::neg(term);
                                    uint32_t id = B->intern(e2.perm, perm_inv(e2.perm, MAXN));
                                    B->get_or_insert(id, tr.v2).push_back(term);
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
                        Expr qe = mono_expr(e.q, esc);
                        long id = -1;
                        for (const VSum& sv : sums) {
                            Expr sumval = SymEngine::mul(sv.second[0], qe);
                            const auto& trs = trans[sv.first];
                            for (size_t t = 0; t < trs.size(); ++t) {
                                const Trans& tr = trs[t];
                                Expr esf = elem_sym_func_q(k, up, e.perm, e.udiff, tr, zidx[index][sv.first][t], esc, yidx);
                                if (esf.is_null()) continue;
                                Expr term = SymEngine::mul(sumval, esf);
                                if (tr.s < 0) term = SymEngine::neg(term);
                                if (id < 0) id = B->intern(e.perm, perm_inv(e.perm, MAXN));
                                B->get_or_insert((uint32_t)id, tr.v2).push_back(term);
                            }
                        }
                    }
                }
            }

            size_t alive = 0;
            for (auto& vec : B->sums) {
                collapse(vec);
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
            if (vec_basic* c = A->find(sid, id_vmu)) result[A->perms[sid]].push_back((*c)[0]);
    }

    PermDict out;
    out.reserve(result.size());
    for (auto& kv : result) {
        Expr s = SymEngine::add(kv.second);
        if (!is_zero(s)) out.push_back({kv.first, s});
    }
    std::sort(out.begin(), out.end(), [](const std::pair<Perm, Expr>& x, const std::pair<Perm, Expr>& y) { return x.first < y.first; });
    return out;
}

int main(int argc, char** argv) {
    CliArgs args = parse_cli(argc, argv, "schubmult_q_double_core", " [--mixed-var]");
    bool same = true;
    for (const std::string& f : args.flags) {
        if (f == "--mixed-var") {
            same = false;
        } else {
            std::fprintf(stderr, "schubmult_q_double_core: unknown option %s (parabolic, nil-Hecke, --slow and --display-positive are not supported)\n", f.c_str());
            return 2;
        }
    }

    ElemSymCache esc(same);
    PermDict coeff_dict;
    coeff_dict.push_back({args.perms[0], SymEngine::one});
    for (size_t i = 1; i < args.perms.size(); ++i) coeff_dict = schubmult_q_double_fast(coeff_dict, args.perms[i], esc);

    // schubmult_q_double prints sorted by (inv, one-line notation)
    std::stable_sort(coeff_dict.begin(), coeff_dict.end(), [](const std::pair<Perm, Expr>& x, const std::pair<Perm, Expr>& y) {
        int ix = perm_inv(x.first, MAXN), iy = perm_inv(y.first, MAXN);
        return ix != iy ? ix < iy : x.first < y.first;
    });
    for (const auto& kv : coeff_dict) std::printf("%s  %s\n", format_perm(kv.first, args.ascode).c_str(), SymEngine::str(*kv.second).c_str());
    return 0;
}
