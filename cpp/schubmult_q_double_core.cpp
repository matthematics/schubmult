// schubmult_q_double_core: C++ port of schubmult.mult.quantum_double.schubmult_q_double_fast
// (products of quantum double Schubert polynomials; no parabolic / nil-Hecke / positivity yet).
//
// Build:  make build/schubmult_q_double_core     (needs SymEngine C++ headers/lib, e.g. from conda)
// Usage:  ./build/schubmult_q_double_core 3 1 2 - 2 1 3                (z = y, the default)
//         ./build/schubmult_q_double_core --mixed-var 3 1 2 - 2 1 3    (S_u(x; y) * S_v(x; z))
//         ./build/schubmult_q_double_core --code 2 0 - 1 0
// Output lines are "(w1, w2, ...)  polynomial", sorted by (inv(w), w), zero coefficients dropped.

#include "kernel_q_double.h"

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
    ExprDict coeff_dict;
    coeff_dict.push_back({args.perms[0], ex_one()});
    for (size_t i = 1; i < args.perms.size(); ++i) coeff_dict = schubmult_q_double_fast(coeff_dict, args.perms[i], esc);

    // schubmult_q_double prints sorted by (inv, one-line notation)
    std::stable_sort(coeff_dict.begin(), coeff_dict.end(), [](const std::pair<Perm, Expr>& x, const std::pair<Perm, Expr>& y) {
        int ix = perm_inv(x.first, MAXN), iy = perm_inv(y.first, MAXN);
        return ix != iy ? ix < iy : x.first < y.first;
    });
    for (const auto& kv : coeff_dict) std::printf("%s  %s\n", format_perm(kv.first, args.ascode).c_str(), ex_str(kv.second).c_str());
    return 0;
}
