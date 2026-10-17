// Scaling experiment: time lrcalc-rs mult_schubert vs the schubmult native kernel on selected pairs.
// Usage: scaling_pairs < pairs.txt  (lines: "u1 .. un - v1 .. vn"), prints "rs_ns sm_cold_ns sm_warm_ns nterms line"
// Built twice is not needed: this binary links the Rust liblrcalc and includes the schubmult kernel.
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <ctime>
#include <sstream>
#include <string>
#include <vector>

#include <lrcalc/ivector.h>
#include <lrcalc/ivlincomb.h>
#include <lrcalc/schublib.h>

#include "../../cpp/kernel_single.h"

static long ns_between(const timespec& a, const timespec& b) { return (b.tv_sec - a.tv_sec) * 1000000000L + (b.tv_nsec - a.tv_nsec); }

int main() {
    std::string line;
    timespec t0, t1;
    while (std::getline(std::cin, line)) {
        std::vector<int> u, v, *cur = &u;
        std::istringstream ss(line);
        std::string tok;
        while (ss >> tok) {
            if (tok == "-") { cur = &v; continue; }
            if (tok[0] == '#') break;
            cur->push_back(std::stoi(tok));
        }
        if (u.empty() || v.empty()) continue;
        int n = (int)std::max(u.size(), v.size());
        while ((int)u.size() < n) u.push_back((int)u.size() + 1);
        while ((int)v.size() < n) v.push_back((int)v.size() + 1);

        ivector *iu = iv_new(n), *iv = iv_new(n);
        for (int i = 0; i < n; i++) { iv_elem(iu, i) = u[i]; iv_elem(iv, i) = v[i]; }
        clock_gettime(CLOCK_MONOTONIC, &t0);
        ivlincomb* res = mult_schubert(iu, iv, 0);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        long rs_ns = ns_between(t0, t1);
        uint32_t nt = 0;
        ivlc_iter itr;
        for (ivlc_first(res, &itr); ivlc_good(&itr); ivlc_next(&itr)) nt++;
        ivlc_free_all(res);
        iv_free(iu); iv_free(iv);

        Perm pu = perm_from_array(u), pv = perm_from_array(v);
        IntDict pd = {{pu, 1}};
        int amb = 2 * n - 1;
        clock_gettime(CLOCK_MONOTONIC, &t0);
        IntDict r1 = schubmult_single(pd, pv, amb);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        long sm_cold = ns_between(t0, t1);
        MultSetup S = mult_setup(pv);
        clock_gettime(CLOCK_MONOTONIC, &t0);
        IntDict r2 = schubmult_single_with_setup(pd, S, amb);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        long sm_warm = ns_between(t0, t1);
        uint32_t nt2 = 0;
        for (auto& kv : r1) if (kv.second) nt2++;
        printf("%ld %ld %ld %u %u %s\n", rs_ns, sm_cold, sm_warm, nt, nt2, line.c_str());
        fflush(stdout);
    }
    return 0;
}
