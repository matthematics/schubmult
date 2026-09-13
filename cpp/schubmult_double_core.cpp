// schubmult_double_core: C++ port of schubmult.mult.double.schubmult_double
// (products of double Schubert polynomials S_u(x; y) * S_v(x; z), expanded in S_w(x; y)),
// with coefficients as SymEngine expressions.
//
// Build:  make build/schubmult_double_core     (needs SymEngine C++ headers/lib, e.g. from conda)
// Usage:  ./build/schubmult_double_core 3 1 2 - 2 1 3                (z = y, the default, as schubmult_double)
//         ./build/schubmult_double_core --mixed-var 3 1 2 - 2 1 3    (S_u(x; y) * S_v(x; z))
//         ./build/schubmult_double_core --code 2 0 - 1 0
// Output lines are "(w1, w2, ...)  polynomial", sorted by (inv(w), w), zero coefficients dropped.

#include "kernel_double.h"
#include "positivity.h"

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
    ExprDict coeff_dict;
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
