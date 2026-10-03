// scaling_feat.cpp: per-pair features and timings of the two native kernels on the scaling families.
// stdin lines: "<anything> u1 .. un - v1 .. vn  # n=N fam"  (scaling.txt or pairs.txt); stdout:
//   n fam inv_u inv_v pd_e nvp_e nterms t_vpath_ns t_trans_ns hybrid_choice(0/1)
// Both kernels work on the factor e with fewer inversions (v-path: recurse on e; transition: expand e).
#include <cstdio>
#include <ctime>
#include <iostream>
#include <sstream>
#include <string>

#include "../../cpp/kernel_transition.h"

static long ns(const timespec& a, const timespec& b) { return (b.tv_sec - a.tv_sec) * 1000000000L + (b.tv_nsec - a.tv_nsec); }

int main() {
    std::string line;
    while (std::getline(std::cin, line)) {
        auto hash = line.find('#');
        if (hash == std::string::npos) continue;
        std::string tag = line.substr(hash + 1), body = line.substr(0, hash);
        std::istringstream ts(tag);
        std::string ntag, fam;
        ts >> ntag >> fam;
        int n = std::stoi(ntag.substr(2));
        // the perms are the last 2n+1 tokens of body
        std::istringstream bs(body);
        std::vector<std::string> toks;
        for (std::string t; bs >> t;) toks.push_back(t);
        std::vector<int> u, v;
        size_t base = toks.size() - (2 * n + 1);
        for (int i = 0; i < n; ++i) u.push_back(std::stoi(toks[base + i]));
        for (int i = 0; i < n; ++i) v.push_back(std::stoi(toks[base + n + 1 + i]));
        Perm pu = perm_from_array(u), pv = perm_from_array(v);
        int amb = 2 * n - 1;
        int inv_u = perm_inv(pu, amb), inv_v = perm_inv(pv, amb);
        const Perm& e = inv_u < inv_v ? pu : pv;
        const Perm& o = inv_u < inv_v ? pv : pu;
        timespec t0, t1;
        clock_gettime(CLOCK_MONOTONIC, &t0);
        IntDict r1 = schubmult_single({{o, 1}}, e, amb);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        long t_vp = ns(t0, t1);
        clock_gettime(CLOCK_MONOTONIC, &t0);
        IntDict r2 = mult_poly_schubert(trans_polynomial(e, amb), o, amb);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        long t_tr = ns(t0, t1);
        MultSetup S = mult_setup(e);
        long nvp = S.trivial ? 0 : vpath_count(S);
        long pd = pipe_dream_count(e, amb);
        MultSetup S2;
        HybridChoice c = hybrid_choose(pu, pv, amb, S2);
        size_t nt = 0;
        for (auto& kv : r1) nt += kv.second != 0;
        size_t nt2 = 0;
        for (auto& kv : r2) nt2 += kv.second != 0;
        if (nt != nt2) fprintf(stderr, "MISMATCH %s\n", line.c_str());
        printf("%d %s %d %d %ld %ld %zu %ld %ld %d\n", n, fam.c_str(), inv_u, inv_v, pd, nvp, nt, t_vp, t_tr, (int)c.transition);
        fflush(stdout);
    }
}
