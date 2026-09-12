// schub_quantum.h: quantum-step combinatorics shared by the quantum kernels
// (schub_lib.elem_sym_perms_q / double_elem_sym_q, perm_utils.count_bruhat).

#pragma once

#include "schub_common.h"

// e[i] = exponent of q_{i+1}
struct Mono {
    uint8_t e[MAXN];
    bool operator==(const Mono& o) const { return std::memcmp(e, o.e, MAXN) == 0; }
    bool operator<(const Mono& o) const { return std::memcmp(e, o.e, MAXN) < 0; }
};

static Mono mono_one() {
    Mono m;
    std::memset(m.e, 0, MAXN);
    return m;
}

static inline Mono mono_mul(const Mono& a, const Mono& b) {
    Mono r;
    for (int i = 0; i < MAXN; ++i) r.e[i] = (uint8_t)(a.e[i] + b.e[i]);
    return r;
}

// perm_utils.count_bruhat
static inline int count_bruhat(const Perm& a, int i, int j) {
    int ai = a.p[i], aj = a.p[j];
    int up = ai < aj ? 1 : -1;
    for (int q = i + 1; q < j; ++q) {
        int aq = a.p[q];
        if (ai < aq && aq < aj)
            up += 2;
        else if (ai > aq && aq > aj)
            up -= 2;
    }
    return up;
}

struct QUp {
    Perm perm;
    int udiff;
    Mono q;
    bool operator==(const QUp& o) const { return udiff == o.udiff && perm == o.perm && q == o.q; }
};

// Appends (perm, udiff, q-monomial) triples to out (cleared first). Positions j run over all
// of [k, MAXN); those beyond the Python's max(k+1, len) bound never satisfy the Bruhat test.
static void elem_sym_perms_q(const Perm& orig, int p, int k, std::vector<QUp>& out) {
    struct E {
        Perm perm;
        Mono q;
        int last_j;
    };
    static std::vector<E> cur, nxt;
    out.clear();
    out.push_back({orig, 0, mono_one()});
    if (p <= 0) return;
    cur.clear();
    cur.push_back({orig, mono_one(), MAXN - 1});
    for (int pp = 0; pp < p; ++pp) {
        nxt.clear();
        for (const E& e : cur) {
            const Perm& up = e.perm;
            int pos[MAXN], npos = 0;
            for (int i = 0; i < k; ++i)
                if (up.p[i] == orig.p[i]) pos[npos++] = i;
            if (npos == 0) continue;
            for (int j = std::min(MAXN - 1, e.last_j); j >= k; --j) {
                for (int t = 0; t < npos; ++t) {
                    int i = pos[t];
                    int ct = count_bruhat(up, i, j);
                    if (ct != 1 && ct != 2 * (i - j) + 1) continue;
                    Mono nq = e.q;
                    if (ct < 0)
                        for (int idx = i; idx < j; ++idx) ++nq.e[idx];  // q_{i+1} ... q_j
                    Perm np = swapped(up, i, j);
                    nxt.push_back({np, nq, j});
                    out.push_back({np, pp + 1, nq});
                }
            }
        }
        std::swap(cur, nxt);
    }
}

// Nontrivial cycles of a permutation: max element (1-indexed) and a membership bitmask.
struct Cycle {
    int mx;
    uint64_t mask[(MAXN + 63) / 64];
};

static void cycles_of(const Perm& s, std::vector<Cycle>& out) {
    out.clear();
    bool seen[MAXN] = {false};
    for (int i = 0; i < MAXN; ++i) {
        if (seen[i] || s.p[i] == i + 1) continue;
        Cycle c;
        std::memset(c.mask, 0, sizeof c.mask);
        c.mx = 0;
        int j = i;
        while (!seen[j]) {
            seen[j] = true;
            c.mask[j >> 6] |= 1ULL << (j & 63);
            c.mx = std::max(c.mx, j + 1);
            j = s.p[j] - 1;
        }
        out.push_back(c);
    }
}

// Python: for each cycle c2 of ip1*perm2 and each cycle c1 of iu*perm1 with the same maximum,
// no other element of c2 may lie in c1.
static bool double_step_good(const std::vector<Cycle>& c1s, const std::vector<Cycle>& c2s) {
    for (const Cycle& c2 : c2s)
        for (const Cycle& c1 : c1s) {
            if (c1.mx != c2.mx) continue;
            int common = 0;
            for (size_t w = 0; w < sizeof c1.mask / 8; ++w) common += __builtin_popcountll(c1.mask[w] & c2.mask[w]);
            if (common > 1) return false;
        }
    return true;
}

// schub_lib.double_elem_sym_q. keys[] are the distinct (perm1, udiff1, q1) of elem_sym_perms_q(u, p1, k);
// second[kk] lists the admissible (perm2, udiff2, q2) of elem_sym_perms_q(perm1, p2, k), accumulated
// over duplicate perm1 entries exactly as the Python dict does.
static void double_elem_sym_q(const Perm& u, int p1, int p2, int k, std::vector<QUp>& keys, std::vector<std::vector<QUp>>& second) {
    static std::vector<QUp> perms1, perms2;
    static std::vector<Cycle> cyc1, cyc2;
    keys.clear();
    second.clear();
    elem_sym_perms_q(u, p1, k, perms1);
    Perm iu = inverse(u, MAXN);
    for (const QUp& e1 : perms1) {
        elem_sym_perms_q(e1.perm, p2, k, perms2);
        cycles_of(compose(iu, e1.perm, MAXN), cyc1);
        Perm ip1 = inverse(e1.perm, MAXN);
        size_t kk = 0;
        while (kk < keys.size() && !(keys[kk] == e1)) ++kk;
        if (kk == keys.size()) {
            keys.push_back(e1);
            second.emplace_back();
        }
        for (const QUp& e2 : perms2) {
            cycles_of(compose(ip1, e2.perm, MAXN), cyc2);
            if (double_step_good(cyc1, cyc2)) second[kk].push_back(e2);
        }
    }
}
