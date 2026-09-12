// schubmult_double_core: C++ port of schubmult.mult.double.schubmult_double
// (products of double Schubert polynomials S_u(x; y) * S_v(x; z), expanded in S_w(x; y)),
// with coefficients as SymEngine expressions.
//
// Build:  make build/schubmult_double_core     (needs SymEngine C++ headers/lib, e.g. from conda)
// Usage:  ./build/schubmult_double_core 3 1 2 - 2 1 3                (z = y, the default, as schubmult_double)
//         ./build/schubmult_double_core --mixed-var 3 1 2 - 2 1 3    (S_u(x; y) * S_v(x; z))
//         ./build/schubmult_double_core --code 2 0 - 1 0
// Output lines are "(w1, w2, ...)  polynomial", sorted by (inv(w), w), zero coefficients dropped.

#include "schub_symbolic.h"
#include "positivity.h"

// ---------------------------------------------------------------------------
// schubmult_double
// ---------------------------------------------------------------------------

typedef std::vector<std::pair<Perm, Expr>> PermDict;

static PermDict schubmult_double(const PermDict& perm_dict, const Perm& v, int n, ElemSymCache& esc) {
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
    typedef PermTable<vec_basic> Table;
    typedef Table::VSum VSum;
    std::unordered_map<Perm, vec_basic, PermHash> result;

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
        A->sums[A->intern(u, inv_u)].push_back({vp.id_start, vec_basic{kv.second}});

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
                                term = SymEngine::mul(sumval, esf);
                            }
                            if (tr.s < 0) term = SymEngine::neg(term);
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
                    Expr s = SymEngine::add(e.second);
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
    CliArgs args = parse_cli(argc, argv, "schubmult_double_core", " [--mixed-var] [--display-positive] [--optimizer-message]");
    bool same = true, display_positive = false, msg = false;
    for (const std::string& f : args.flags) {
        if (f == "--mixed-var") {
            same = false;
        } else if (f == "--display-positive") {
            display_positive = true;
        } else if (f == "--optimizer-message") {
            msg = true;
        } else {
            std::fprintf(stderr, "schubmult_double_core: unknown option %s\n", f.c_str());
            return 2;
        }
    }
    if (display_positive && msg && same) std::fprintf(stderr, "schubmult_double_core: --optimizer-message has no effect without --mixed-var\n");

    ElemSymCache esc(same);
    // Order matters here (u carries y, v carries z), so multiply in the order given.
    PermDict coeff_dict;
    coeff_dict.push_back({args.perms[0], SymEngine::one});
    for (size_t i = 1; i < args.perms.size(); ++i) coeff_dict = schubmult_double(coeff_dict, args.perms[i], args.n, esc);

    // Mixed variables: every coefficient goes straight to the MILP (the posify shortcut formulas are
    // not ported yet). Same variables: sv_posify rewrites in the roots y_{j+1} - y_j.
    if (display_positive) {
        for (auto& kv : coeff_dict) kv.second = same ? positivity::sv_posify(kv.second, esc.Y) : positivity::compute_positive_rep(kv.second, esc.Y, esc.Z, msg);
        coeff_dict.erase(std::remove_if(coeff_dict.begin(), coeff_dict.end(), [](const std::pair<Perm, Expr>& kv) { return is_zero(kv.second); }), coeff_dict.end());
    }

    // schubmult_double prints sorted by (inv, one-line notation), right-aligned to the widest permutation
    std::stable_sort(coeff_dict.begin(), coeff_dict.end(), [](const std::pair<Perm, Expr>& x, const std::pair<Perm, Expr>& y) {
        int ix = perm_inv(x.first, MAXN), iy = perm_inv(y.first, MAXN);
        return ix != iy ? ix < iy : x.first < y.first;
    });
    size_t width = 0;
    for (const auto& kv : coeff_dict) width = std::max(width, format_perm(kv.first, args.ascode).size());
    for (const auto& kv : coeff_dict) std::printf("%*s  %s\n", (int)width, format_perm(kv.first, args.ascode).c_str(), SymEngine::str(*kv.second).c_str());
    return 0;
}
