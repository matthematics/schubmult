// schub_common.h: permutation utilities, v-path precomputation, elem_sym_perms,
// and the intern table shared by schubmult_core (single), schubmult_double_core and
// schubmult_q_core.
//
// Every permutation is one fixed-size byte array of length MAXN padded with the identity,
// so equal permutations are byte-identical; rebuild with -DMAXN=64 for larger inputs
// (u,v in S_a,S_b need a+b-1 <= MAXN and 2b-1 <= MAXN).
#pragma once

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

#ifndef MAXN
#define MAXN 32
#endif
static_assert(MAXN % 8 == 0 && MAXN <= 248, "MAXN must be a multiple of 8 and at most 248");

// ---------------------------------------------------------------------------
// Permutation type
// ---------------------------------------------------------------------------

struct Perm {
    uint8_t p[MAXN];
    bool operator==(const Perm& o) const { return std::memcmp(p, o.p, MAXN) == 0; }
    bool operator!=(const Perm& o) const { return !(*this == o); }
    bool operator<(const Perm& o) const { return std::memcmp(p, o.p, MAXN) < 0; }
};

struct PermHash {
    size_t operator()(const Perm& x) const {
        uint64_t h = 0x9E3779B97F4A7C15ULL;
        for (int i = 0; i < MAXN; i += 8) {
            uint64_t w;
            std::memcpy(&w, x.p + i, 8);
            h ^= w;
            h *= 0xFF51AFD7ED558CCDULL;
            h ^= h >> 33;
        }
        return (size_t)h;
    }
};

[[noreturn]] static void die(const char* msg) {
    std::fprintf(stderr, "schubmult: %s\n", msg);
    std::exit(1);
}

static Perm identity_perm() {
    Perm r;
    for (int i = 0; i < MAXN; ++i) r.p[i] = (uint8_t)(i + 1);
    return r;
}

static int perm_inv(const Perm& a, int n) {
    int c = 0;
    for (int i = 0; i < n; ++i) {
        int ai = a.p[i];
        for (int j = i + 1; j < n; ++j) c += ai > a.p[j];
    }
    return c;
}

// (a * b)[i] = a[b[i] - 1]  (matches Permutation.__mul__)
static Perm compose(const Perm& a, const Perm& b, int n) {
    Perm r = identity_perm();
    for (int i = 0; i < n; ++i) r.p[i] = a.p[b.p[i] - 1];
    return r;
}

static Perm inverse(const Perm& a, int n) {
    Perm r = identity_perm();
    for (int i = 0; i < n; ++i) r.p[a.p[i] - 1] = (uint8_t)(i + 1);
    return r;
}

static std::vector<int> lehmer_code(const Perm& a, int n) {
    std::vector<int> cd(n, 0);
    for (int i = 0; i < n; ++i) {
        int ai = a.p[i], c = 0;
        for (int j = i + 1; j < n; ++j) c += a.p[j] < ai;
        cd[i] = c;
    }
    return cd;
}

// Permutation._cached_theta, with trailing zeros removed.
static std::vector<int> theta(std::vector<int> cd) {
    int L = (int)cd.size();
    for (int i = L - 1; i > 0; --i)
        for (int j = i - 1; j >= 0; --j)
            if (cd[j] < cd[i]) cd[i] += 1;
    std::sort(cd.begin(), cd.end(), [](int x, int y) { return x > y; });
    while (!cd.empty() && cd.back() == 0) cd.pop_back();
    return cd;
}

// Permutation._cached_medium_theta (used by schubmult_q_fast), with trailing zeros removed.
static std::vector<int> medium_theta(std::vector<int> cd) {
    int L = (int)cd.size();
    bool found_one = true;
    while (found_one) {
        found_one = false;
        for (int i = 0; i + 1 < L; ++i) {
            if (cd[i] < cd[i + 1]) {
                int a = cd[i], b = cd[i + 1];
                cd[i] = b + 1;
                cd[i + 1] = a;
                found_one = true;
                break;
            }
            if (cd[i] == cd[i + 1] && cd[i] != 0 && i > 0 && cd[i - 1] <= cd[i] + 1) {
                cd[i] += 1;
                found_one = true;
                break;
            }
        }
    }
    while (!cd.empty() && cd.back() == 0) cd.pop_back();
    return cd;
}

static Perm uncode(const std::vector<int>& cd_in) {
    if (cd_in.empty()) return identity_perm();
    int L = (int)cd_in.size(), maxreq = 0;
    for (int i = 0; i < L; ++i) maxreq = std::max(maxreq, cd_in[i] + i);
    std::vector<int> cd(cd_in);
    if (maxreq > L) cd.resize(maxreq, 0);
    int len = (int)cd.size() + 1;
    if (len > MAXN) die("permutation too long for MAXN; rebuild with a larger -DMAXN");
    std::vector<int> full(len);
    for (int i = 0; i < len; ++i) full[i] = i + 1;
    Perm r = identity_perm();
    for (size_t i = 0; i < cd.size(); ++i) {
        r.p[i] = (uint8_t)full[cd[i]];
        full.erase(full.begin() + cd[i]);
    }
    r.p[cd.size()] = (uint8_t)full[0];
    return r;
}

static inline Perm swapped(const Perm& a, int i, int j) {
    Perm r = a;
    std::swap(r.p[i], r.p[j]);
    return r;
}

// i < j. a[i] < a[j] and no value strictly between them in positions i+1..j-1.
static inline bool bruhat_ascent(const Perm& a, int i, int j) {
    int ai = a.p[i], aj = a.p[j];
    if (ai > aj) return false;
    for (int q = i + 1; q < j; ++q) {
        int aq = a.p[q];
        if (ai < aq && aq < aj) return false;
    }
    return true;
}

// i < j. a[i] > a[j] and no value strictly between them in positions i+1..j-1.
static inline bool bruhat_descent(const Perm& a, int i, int j) {
    int ai = a.p[i], aj = a.p[j];
    if (ai < aj) return false;
    for (int q = i + 1; q < j; ++q) {
        int aq = a.p[q];
        if (aq > aj && ai > aq) return false;
    }
    return true;
}

// ---------------------------------------------------------------------------
// v-side precomputation (compute_vpathdicts / kdown_perms)
// ---------------------------------------------------------------------------

struct Trip {
    Perm perm;
    int pp;
    int s;
    bool operator<(const Trip& o) const {
        if (perm != o.perm) return perm < o.perm;
        if (pp != o.pp) return pp < o.pp;
        return s < o.s;
    }
    bool operator==(const Trip& o) const { return perm == o.perm && pp == o.pp && s == o.s; }
};

static std::vector<Trip> kdown_perms(const Perm& perm, const Perm& monoperm, int p, int k) {
    const int n = MAXN;
    int inv_m = perm_inv(monoperm, n), inv_p = perm_inv(perm, n);
    std::vector<Trip> full;
    if (perm_inv(compose(perm, monoperm, n), n) == inv_m - inv_p) full.push_back({perm, 0, 1});
    int a2 = k - 1;
    if (a2 >= n) return full;
    std::vector<std::pair<Perm, int>> down, down2;
    down.push_back({perm, 1});
    for (int pp = 1; pp <= p; ++pp) {
        down2.clear();
        int g_inv = inv_m - inv_p + pp;
        for (const auto& e : down) {
            const Perm& perm2 = e.first;
            int s = e.second, s2 = -s;
            for (int b = 0; b < k - 1; ++b) {
                if (perm2.p[b] != perm.p[b] || !bruhat_descent(perm2, b, a2)) continue;
                Perm np = swapped(perm2, b, a2);
                down2.push_back({np, s2});
                if (perm_inv(compose(np, monoperm, n), n) == g_inv) full.push_back({np, pp, s2});
            }
            for (int b = k; b < n; ++b) {
                if (perm2.p[b] != perm.p[b] || !bruhat_descent(perm2, a2, b)) continue;
                Perm np = swapped(perm2, a2, b);
                down2.push_back({np, s});
                if (perm_inv(compose(np, monoperm, n), n) == g_inv) full.push_back({np, pp, s});
            }
        }
        // Duplicate (perm, s) chains generate identical descendants; the result is a set anyway.
        std::sort(down2.begin(), down2.end());
        down2.erase(std::unique(down2.begin(), down2.end()), down2.end());
        std::swap(down, down2);
    }
    return full;
}

struct Trans {
    uint32_t v2;
    int vdiff;
    int s;
    bool operator<(const Trans& o) const {
        if (v2 != o.v2) return v2 < o.v2;
        if (vdiff != o.vdiff) return vdiff < o.vdiff;
        return s < o.s;
    }
    bool operator==(const Trans& o) const { return v2 == o.v2 && vdiff == o.vdiff && s == o.s; }
};

struct VPaths {
    // level[i] = v-perms reachable after i steps; level[0] = {id}, level[thL] = {vmu}
    std::vector<std::vector<Perm>> level;
    // trans[i][vp_id] = list of (v2_id in level[i+1], vdiff, s)
    std::vector<std::vector<std::vector<Trans>>> trans;
    std::vector<int> mx_th;
    uint32_t id_start;  // id of identity in level[0], or UINT32_MAX
};

static VPaths compute_vpathdicts(const std::vector<int>& th, const Perm& vmu) {
    int thL = (int)th.size();
    std::vector<std::unordered_map<Perm, std::vector<Trip>, PermHash>> vpd(thL);
    vpd[thL - 1][vmu] = {};
    for (int i = thL - 1; i >= 0; --i) {
        std::vector<int> prefix(th.begin(), th.begin() + i);
        Perm monoperm = inverse(uncode(prefix), MAXN);
        int k = i + 1;
        std::vector<Perm> keys;
        keys.reserve(vpd[i].size());
        for (auto& kv : vpd[i]) keys.push_back(kv.first);
        for (const Perm& last : keys) {
            std::vector<Trip> trips = kdown_perms(last, monoperm, th[i], k);
            std::sort(trips.begin(), trips.end());
            trips.erase(std::unique(trips.begin(), trips.end()), trips.end());
            if (i > 0)
                for (const Trip& t : trips) vpd[i - 1].try_emplace(t.perm);
            vpd[i][last] = std::move(trips);
        }
    }

    VPaths vp;
    vp.level.resize(thL + 1);
    for (int i = 0; i < thL; ++i) {
        for (auto& kv : vpd[i]) vp.level[i + 1].push_back(kv.first);
        std::sort(vp.level[i + 1].begin(), vp.level[i + 1].end());
    }
    for (auto& kv : vpd[0])
        for (const Trip& t : kv.second) vp.level[0].push_back(t.perm);
    std::sort(vp.level[0].begin(), vp.level[0].end());
    vp.level[0].erase(std::unique(vp.level[0].begin(), vp.level[0].end()), vp.level[0].end());

    std::vector<std::unordered_map<Perm, uint32_t, PermHash>> ids(thL + 1);
    for (int i = 0; i <= thL; ++i)
        for (uint32_t j = 0; j < vp.level[i].size(); ++j) ids[i][vp.level[i][j]] = j;

    vp.trans.resize(thL);
    vp.mx_th.assign(thL, 0);
    for (int i = 0; i < thL; ++i) {
        vp.trans[i].resize(vp.level[i].size());
        for (auto& kv : vpd[i]) {
            uint32_t v2 = ids[i + 1][kv.first];
            for (const Trip& t : kv.second) {
                uint32_t vpid = ids[i][t.perm];
                vp.trans[i][vpid].push_back({v2, t.pp, t.s});
                vp.mx_th[i] = std::max(vp.mx_th[i], th[i] - t.pp);
            }
        }
        for (auto& lst : vp.trans[i]) {
            std::sort(lst.begin(), lst.end());
            lst.erase(std::unique(lst.begin(), lst.end()), lst.end());
        }
    }
    auto it = ids[0].find(identity_perm());
    vp.id_start = it == ids[0].end() ? UINT32_MAX : it->second;
    return vp;
}

// ---------------------------------------------------------------------------
// u-side: elem_sym_perms
// ---------------------------------------------------------------------------

struct UpEntry {
    Perm perm;
    int last;
};

// Appends (perm, udiff) pairs to out (cleared first). n bounds the positions considered.
[[maybe_unused]] static void elem_sym_perms(const Perm& orig, int p, int k, int n, std::vector<std::pair<Perm, int>>& out) {
    static std::vector<UpEntry> cur, nxt;
    out.clear();
    out.push_back({orig, 0});
    if (p <= 0) return;
    cur.clear();
    cur.push_back({orig, 1 << 30});
    int kk = std::min(k, n);
    for (int pp = 0; pp < p; ++pp) {
        nxt.clear();
        for (const UpEntry& e : cur) {
            const Perm& up = e.perm;
            int last = e.last;
            int pos[MAXN], npos = 0;
            for (int i = 0; i < kk; ++i)
                if (up.p[i] < last) pos[npos++] = i;
            if (npos == 0) continue;
            for (int j = k; j < n; ++j) {
                int uj = up.p[j];
                if (uj >= last) continue;
                for (int t = 0; t < npos; ++t) {
                    int i = pos[t];
                    if (!bruhat_ascent(up, i, j)) continue;
                    Perm np = swapped(up, i, j);
                    nxt.push_back({np, uj});
                    out.push_back({np, pp + 1});
                }
            }
        }
        std::swap(cur, nxt);
    }
}

// ---------------------------------------------------------------------------
// Intern table: every distinct permutation of one step is stored exactly once and
// referred to by a dense 4-byte id. Open addressing, linear probing. Carries the
// inversion number and the v-path sums (pairs (v_id, T)) of every permutation.
// ---------------------------------------------------------------------------

template <typename T>
struct PermTable {
    typedef std::pair<uint32_t, T> VSum;
    std::vector<Perm> perms;
    std::vector<int> inv;
    std::vector<std::vector<VSum>> sums;
    std::vector<uint32_t> slots;  // id + 1, 0 = empty
    uint32_t count = 0;

    void reset() {
        perms.clear();
        inv.clear();
        sums.clear();
        count = 0;
        if (slots.size() < 64)
            slots.assign(64, 0);
        else
            std::fill(slots.begin(), slots.end(), 0);
    }

    void grow() {
        std::vector<uint32_t> ns(slots.size() * 2, 0);
        size_t mask = ns.size() - 1;
        for (uint32_t id = 0; id < count; ++id) {
            size_t i = PermHash()(perms[id]) & mask;
            while (ns[i]) i = (i + 1) & mask;
            ns[i] = id + 1;
        }
        slots.swap(ns);
    }

    // Returns the id of p, inserting it (with inversion number inv_p) if new.
    uint32_t intern(const Perm& p, int inv_p) {
        size_t mask = slots.size() - 1;
        size_t i = PermHash()(p) & mask;
        for (;;) {
            uint32_t s = slots[i];
            if (s == 0) break;
            if (perms[s - 1] == p) return s - 1;
            i = (i + 1) & mask;
        }
        uint32_t id = count++;
        slots[i] = id + 1;
        perms.push_back(p);
        inv.push_back(inv_p);
        sums.emplace_back();
        if ((size_t)count * 2 > slots.size()) grow();
        return id;
    }

    // Pointer to the T stored for v_id (a few entries each, so a linear scan wins), or nullptr.
    T* find(uint32_t id, uint32_t v) {
        for (VSum& e : sums[id])
            if (e.first == v) return &e.second;
        return nullptr;
    }
    T& get_or_insert(uint32_t id, uint32_t v) {
        if (T* t = find(id, v)) return *t;
        sums[id].push_back({v, T()});
        return sums[id].back().second;
    }
};

#ifdef STATS
static double now_s() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
#endif

// ---------------------------------------------------------------------------
// Setup shared by the multiplication kernels: theta, mu, vmu, v-paths.
// ---------------------------------------------------------------------------

struct MultSetup {
    std::vector<int> th;
    Perm mu, vmu;
    int inv_mu = 0, inv_vmu = 0;
    VPaths vp;
    bool trivial = false;  // v is the identity
};

static MultSetup mult_setup(const Perm& v, bool medium = false) {
    MultSetup s;
    std::vector<int> cd = lehmer_code(inverse(v, MAXN), MAXN);
    s.th = medium ? medium_theta(cd) : theta(cd);
    if (s.th.empty()) {
        s.trivial = true;
        return s;
    }
    s.mu = uncode(s.th);
    s.vmu = compose(v, s.mu, MAXN);
    s.inv_vmu = perm_inv(s.vmu, MAXN);
    s.inv_mu = perm_inv(s.mu, MAXN);
#ifdef STATS
    double t_pre = now_s();
#endif
    s.vp = compute_vpathdicts(s.th, s.vmu);
#ifdef STATS
    std::fprintf(stderr, "th=[");
    for (int t : s.th) std::fprintf(stderr, "%d ", t);
    std::fprintf(stderr, "] vpathdicts %.3fs; level sizes:", now_s() - t_pre);
    for (auto& L : s.vp.level) std::fprintf(stderr, " %zu", L.size());
    std::fprintf(stderr, "\n");
#endif
    return s;
}

// ---------------------------------------------------------------------------
// CLI helpers
// ---------------------------------------------------------------------------

static std::string format_perm(const Perm& a, bool ascode) {
    std::string s;
    if (ascode) {
        std::vector<int> cd = lehmer_code(a, MAXN);
        while (!cd.empty() && cd.back() == 0) cd.pop_back();
        s = "[";
        for (size_t i = 0; i < cd.size(); ++i) {
            if (i) s += ", ";
            s += std::to_string(cd[i]);
        }
        s += "]";
    } else {
        int L = MAXN;
        while (L > 0 && a.p[L - 1] == L) --L;
        s = "(";
        for (int i = 0; i < L; ++i) {
            if (i) s += ", ";
            s += std::to_string((int)a.p[i]);
        }
        s += ")";
    }
    return s;
}

static Perm perm_from_array(const std::vector<int>& arr) {
    int L = (int)arr.size();
    if (L > MAXN) die("permutation too long for MAXN; rebuild with a larger -DMAXN");
    std::vector<char> seen(L + 1, 0);
    for (int x : arr) {
        if (x < 1 || x > L || seen[x]) die("input is not a permutation of 1..n");
        seen[x] = 1;
    }
    Perm r = identity_perm();
    for (int i = 0; i < L; ++i) r.p[i] = (uint8_t)arr[i];
    return r;
}

struct CliArgs {
    bool ascode = false;
    std::vector<std::string> flags;  // other --options, for the caller
    std::vector<Perm> perms;
    int n = 2;  // every permutation in the product lies in S_n (S_a * S_b lives in S_{a+b-1})
};

static CliArgs parse_cli(int argc, char** argv, const char* prog, const char* extra_help) {
    CliArgs a;
    std::vector<std::vector<int>> groups(1);
    auto usage = [&]() {
        std::fprintf(stderr,
                     "usage: %s [--code]%s p1 p2 ... - q1 q2 ... [- ...]\n"
                     "  Permutations in one-line notation, or Lehmer codes with --code.\n",
                     prog, extra_help);
        std::exit(2);
    };
    for (int i = 1; i < argc; ++i) {
        std::string s = argv[i];
        if (s == "--code") {
            a.ascode = true;
        } else if (s == "-") {
            groups.emplace_back();
        } else if (s == "-h" || s == "--help") {
            usage();
        } else if (s.size() > 2 && s[0] == '-' && s[1] == '-') {
            a.flags.push_back(s);
        } else {
            char* end = nullptr;
            long x = std::strtol(argv[i], &end, 10);
            if (!end || *end != '\0' || x < 0) {
                std::fprintf(stderr, "%s: bad argument '%s'\n", prog, argv[i]);
                usage();
            }
            groups.back().push_back((int)x);
        }
    }
    if (groups[0].empty()) usage();
    int n = 1;
    for (const auto& g : groups) {
        int len = (int)g.size();
        if (a.ascode) {
            len = g.empty() ? 0 : (int)g.size() + 1;
            for (size_t i = 0; i < g.size(); ++i) len = std::max(len, g[i] + (int)i + 1);
        }
        n += std::max(len, 1) - 1;
        a.perms.push_back(a.ascode ? uncode(g) : perm_from_array(g));
    }
    a.n = std::max(n, 2);
    if (a.n > MAXN) die("product needs S_n with n > MAXN; rebuild with a larger -DMAXN");
    return a;
}
