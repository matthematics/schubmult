// schubmult_q_core: C++ port of schubmult.mult.quantum.schubmult_q_fast
// (products of quantum Schubert polynomials; no parabolic support).
//
// Build:  make build/schubmult_q_core
// Usage:  ./build/schubmult_q_core 3 1 2 - 2 1 3
//         ./build/schubmult_q_core --code 2 0 - 1 0
// Output lines are "(w1, w2, ...)  polynomial in q_i", sorted by (inv(w), w), zeros dropped.
//
// Coefficients are polynomials in q_1, q_2, ... with int coefficients. Quantum steps can
// shrink a permutation, so unlike the classical kernel every position < MAXN is in play.

#include "kernel_q.h"

int main(int argc, char** argv) {
    CliArgs args = parse_cli(argc, argv, "schubmult_q_core", "");
    for (const std::string& f : args.flags) {
        std::fprintf(stderr, "schubmult_q_core: unknown option %s (parabolic and --slow are not supported)\n", f.c_str());
        return 2;
    }

    // schubmult_q multiplies in the order given.
    QDict coeff_dict;
    QPoly one;
    one.t.push_back({mono_one(), 1});
    coeff_dict.push_back({args.perms[0], one});
    for (size_t i = 1; i < args.perms.size(); ++i) coeff_dict = schubmult_q_fast(coeff_dict, args.perms[i]);

    // schubmult_q prints sorted by (inv, one-line notation)
    std::stable_sort(coeff_dict.begin(), coeff_dict.end(), [](const std::pair<Perm, QPoly>& x, const std::pair<Perm, QPoly>& y) {
        int ix = perm_inv(x.first, MAXN), iy = perm_inv(y.first, MAXN);
        return ix != iy ? ix < iy : x.first < y.first;
    });
    for (const auto& kv : coeff_dict) std::printf("%s  %s\n", format_perm(kv.first, args.ascode).c_str(), qpoly_str(kv.second).c_str());
    return 0;
}
