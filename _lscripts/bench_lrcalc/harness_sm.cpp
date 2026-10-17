// Native harness for the schubmult C++ single kernel (same record format as harness.c).
// Usage: harness_sm n start_u end_u outfile [mode]
//   mode 0 (default): time mult_setup(v) + kernel, i.e. one cold call per pair (what the CLI does)
//   mode 1: time the kernel only; the v-path tables of v are built outside the timed region
//           (what a caller multiplying many things by the same v, or the cached Python path, pays)
//   mode 2: swap operands so the recursion runs on the factor with fewer inversions (like lrcalc), cold
//   mode 3: cost-model hybrid (kernel_transition.h schubmult_single_hybrid), cold
//   mode 4: pure transition + Monk, expanding the factor with fewer inversions, cold
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <vector>

#include "../../cpp/kernel_transition.h"

static uint64_t fnv_perm(const Perm& w, int n) {
    int len = n;
    while (len > 0 && w.p[len - 1] == len) len--;
    uint64_t h = 1469598103934665603ULL;
    for (int i = 0; i < len; i++) {
        h ^= (uint64_t)w.p[i];
        h *= 1099511628211ULL;
    }
    return h;
}

static int next_perm(std::vector<int>& a) {
    int n = (int)a.size(), i = n - 2;
    while (i >= 0 && a[i] > a[i + 1]) i--;
    if (i < 0) return 0;
    int j = n - 1;
    while (a[j] < a[i]) j--;
    std::swap(a[i], a[j]);
    for (int l = i + 1, r = n - 1; l < r; l++, r--) std::swap(a[l], a[r]);
    return 1;
}

struct rec_t { uint32_t ui, vi, time_ns, nterms; uint64_t hash; };

int main(int argc, char** argv) {
    if (argc != 5 && argc != 6) { fprintf(stderr, "usage: %s n start_u end_u outfile [mode]\n", argv[0]); return 2; }
    int n = atoi(argv[1]);
    long start_u = atol(argv[2]), end_u = atol(argv[3]);
    int mode = argc == 6 ? atoi(argv[5]) : 0;
    FILE* out = fopen(argv[4], "wb");
    if (!out) { perror("fopen"); return 1; }

    std::vector<Perm> perms;
    std::vector<int> cur(n);
    for (int i = 0; i < n; i++) cur[i] = i + 1;
    do { perms.push_back(perm_from_array(cur)); } while (next_perm(cur));
    long total = (long)perms.size();
    // the product of two elements of S_n lies in S_{2n-1}
    int amb = 2 * n - 1;
    std::vector<int> invs(total);
    for (long i = 0; i < total; i++) invs[i] = perm_inv(perms[i], n);

    std::vector<rec_t> buf; buf.reserve(4096);
    timespec t0, t1;
    for (long ui = start_u; ui < end_u && ui < total; ui++) {
        for (long vi = ui; vi < total; vi++) {
            long a = ui, b = vi;
            if ((mode == 2 || mode == 4) && invs[ui] < invs[vi]) { a = vi; b = ui; }  // recurse on / expand the factor with fewer inversions
            IntDict pd = {{perms[a], 1}};
            IntDict res;
            if (mode == 1) {
                MultSetup S = mult_setup(perms[b]);
                clock_gettime(CLOCK_MONOTONIC, &t0);
                res = schubmult_single_with_setup(pd, S, amb);
                clock_gettime(CLOCK_MONOTONIC, &t1);
            } else if (mode == 3) {
                clock_gettime(CLOCK_MONOTONIC, &t0);
                res = schubmult_single_hybrid(perms[ui], perms[vi], amb);
                clock_gettime(CLOCK_MONOTONIC, &t1);
            } else if (mode == 4) {
                clock_gettime(CLOCK_MONOTONIC, &t0);
                res = schubmult_transition(pd, perms[b], amb);
                clock_gettime(CLOCK_MONOTONIC, &t1);
            } else {
                clock_gettime(CLOCK_MONOTONIC, &t0);
                res = schubmult_single(pd, perms[b], amb);
                clock_gettime(CLOCK_MONOTONIC, &t1);
            }
            uint64_t h = 0; uint32_t nt = 0;
            for (const auto& kv : res) {
                if (kv.second == 0) continue;
                h += fnv_perm(kv.first, amb) * (uint64_t)(int64_t)kv.second;
                nt++;
            }
            long ns = (t1.tv_sec - t0.tv_sec) * 1000000000L + (t1.tv_nsec - t0.tv_nsec);
            buf.push_back({(uint32_t)ui, (uint32_t)vi, (uint32_t)(ns > 0xFFFFFFFFL ? 0xFFFFFFFFL : ns), nt, h});
            if (buf.size() == 4096) { fwrite(buf.data(), sizeof(rec_t), buf.size(), out); buf.clear(); }
        }
    }
    if (!buf.empty()) fwrite(buf.data(), sizeof(rec_t), buf.size(), out);
    fclose(out);
    return 0;
}
