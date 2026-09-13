// schubmult_core: standalone C++ port of schubmult.mult.single.schubmult_py
// (products of ordinary Schubert polynomials via theta codes / v-path dicts).
//
// Build:  make            (binaries go to build/; or: g++ -O3 -march=native -std=c++17 -o build/schubmult_core schubmult_core.cpp)
// Usage:  ./build/schubmult_core 3 1 2 - 2 1 3
//         ./build/schubmult_core --code 2 0 - 1 0
// Output lines are "coeff  (w1, w2, ...)" like schubmult_py (or "[c1, c2]" with --code).

#include "kernel_single.h"

int main(int argc, char** argv) {
    CliArgs args = parse_cli(argc, argv, "schubmult_core", "");
    for (const std::string& f : args.flags) {
        std::fprintf(stderr, "schubmult_core: unknown option %s\n", f.c_str());
        return 2;
    }

    // Same ordering as schubmult_py: largest inv(v*mu) first (stable).
    std::vector<std::pair<int, Perm>> keyed;
    for (const Perm& x : args.perms) {
        std::vector<int> th = theta(lehmer_code(inverse(x, MAXN), MAXN));
        int key = 0;
        for (int t : th) key += t;
        key -= perm_inv(x, MAXN);
        keyed.push_back({key, x});
    }
    std::stable_sort(keyed.begin(), keyed.end(), [](const std::pair<int, Perm>& a, const std::pair<int, Perm>& b) { return a.first > b.first; });

    IntDict coeff_dict;
    coeff_dict.push_back({keyed[0].second, 1});
    for (size_t i = 1; i < keyed.size(); ++i) coeff_dict = schubmult_single(coeff_dict, keyed[i].second, args.n);

    for (const auto& kv : coeff_dict)
        if (kv.second != 0) std::printf("%d  %s\n", kv.second, format_perm(kv.first, args.ascode).c_str());
    return 0;
}
