// positivity.h: C++ port of schubmult.mult.positivity.compute_positive_rep.
//
// Given a coefficient val in y's and z's, finds nonnegative integers a_b with
//     val = sum_b a_b * prod_{(k,i) in b} (y_k - z_i)
// over a combinatorially generated family of root products b, by solving an integer
// feasibility problem with the `cbc` binary (as PuLP's PULP_CBC_CMD does). This is the one
// place where expansion is required: the MILP works on the monomial coefficient vector.

#pragma once

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <unistd.h>
#include <unordered_map>
#include <vector>

#include <symengine/add.h>
#include <symengine/basic.h>
#include <symengine/constants.h>
#include <symengine/integer.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/subs.h>
#include <symengine/symbol.h>
#include <symengine/visitor.h>

namespace positivity {

using SymEngine::Basic;
using SymEngine::RCP;
typedef RCP<const Basic> Expr;

typedef std::vector<int> Exps;  // exponent vector over the variable list (y's first, then z's)

struct ExpsHash {
    size_t operator()(const Exps& v) const {
        uint64_t h = 0x9E3779B97F4A7C15ULL;
        for (int x : v) {
            h ^= (uint64_t)(uint32_t)x;
            h *= 0xFF51AFD7ED558CCDULL;
            h ^= h >> 33;
        }
        return (size_t)h;
    }
};
typedef std::unordered_map<Exps, long, ExpsHash> SPoly;

[[noreturn]] static void fail(const std::string& msg) {
    std::fprintf(stderr, "positivity: %s\n", msg.c_str());
    std::exit(1);
}

static long number_to_long(const Basic& n) {
    if (!SymEngine::is_a<SymEngine::Integer>(n)) fail("non-integer coefficient " + n.__str__());
    return SymEngine::down_cast<const SymEngine::Integer&>(n).as_int();
}

// Adds coef * mono to P, where mono is a Symbol, Pow(Symbol, Integer), or Mul of those.
static void add_monomial(SPoly& P, const Basic& mono, long coef, const std::map<std::string, int>& varidx, int nvars) {
    Exps e(nvars, 0);
    auto put = [&](const Basic& base, long exp) {
        if (!SymEngine::is_a<SymEngine::Symbol>(base)) fail("unexpected factor " + base.__str__());
        auto it = varidx.find(SymEngine::down_cast<const SymEngine::Symbol&>(base).get_name());
        if (it == varidx.end()) fail("unknown variable " + base.__str__());
        e[it->second] += (int)exp;
    };
    if (SymEngine::is_a<SymEngine::Symbol>(mono)) {
        put(mono, 1);
    } else if (SymEngine::is_a<SymEngine::Pow>(mono)) {
        const auto& pw = SymEngine::down_cast<const SymEngine::Pow&>(mono);
        put(*pw.get_base(), number_to_long(*pw.get_exp()));
    } else if (SymEngine::is_a<SymEngine::Mul>(mono)) {
        const auto& m = SymEngine::down_cast<const SymEngine::Mul&>(mono);
        coef *= number_to_long(*m.get_coef());
        for (const auto& kv : m.get_dict()) put(*kv.first, number_to_long(*kv.second));
    } else if (SymEngine::is_a_Number(mono)) {
        coef *= number_to_long(mono);
    } else {
        fail("unexpected term " + mono.__str__());
    }
    P[e] += coef;
}

static SPoly to_spoly(const Expr& expanded, const std::map<std::string, int>& varidx, int nvars) {
    SPoly P;
    if (SymEngine::is_a<SymEngine::Add>(*expanded)) {
        const auto& a = SymEngine::down_cast<const SymEngine::Add&>(*expanded);
        long c0 = number_to_long(*a.get_coef());
        if (c0 != 0) P[Exps(nvars, 0)] += c0;
        for (const auto& kv : a.get_dict()) add_monomial(P, *kv.first, number_to_long(*kv.second), varidx, nvars);
    } else {
        add_monomial(P, *expanded, 1, varidx, nvars);
    }
    for (auto it = P.begin(); it != P.end();)
        it = it->second == 0 ? P.erase(it) : std::next(it);
    return P;
}

// P * (y_k - z_i), with k, i positions in the variable list
static SPoly mul_root(const SPoly& P, int kpos, int ipos) {
    SPoly R;
    for (const auto& kv : P) {
        Exps e = kv.first;
        ++e[kpos];
        R[e] += kv.second;
        --e[kpos];
        ++e[ipos];
        R[e] -= kv.second;
    }
    for (auto it = R.begin(); it != R.end();)
        it = it->second == 0 ? R.erase(it) : std::next(it);
    return R;
}

// Drop monomials divisible by z_1 (the Python's subs({z_1: 0})).
static SPoly kill_z1(const SPoly& P, int z1pos) {
    if (z1pos < 0) return P;
    SPoly R;
    for (const auto& kv : P)
        if (kv.first[z1pos] == 0) R.insert(kv);
    return R;
}

static int index_of(const std::string& name, char prefix) {
    if (name.size() < 3 || name[0] != prefix || name[1] != '_') return -1;
    return std::atoi(name.c_str() + 2);
}

struct BaseVector {
    std::vector<std::pair<int, int>> factors;  // (y position, z position) in the variable list
    SPoly vec;                                 // expanded, z_1 -> 0
};

// The cbc binary: $SCHUBMULT_CBC, else PuLP's bundled copy in the active conda env, else PATH.
static std::string cbc_path() {
    if (const char* e = std::getenv("SCHUBMULT_CBC")) return e;
    if (const char* prefix = std::getenv("CONDA_PREFIX")) {
        std::string glob_dir = std::string(prefix) + "/lib";
        for (const char* py : {"python3.13", "python3.12", "python3.11", "python3.10"}) {
            std::string p = glob_dir + "/" + py + "/site-packages/pulp/solverdir/cbc/linux/i64/cbc";
            if (access(p.c_str(), X_OK) == 0) return p;
        }
    }
    return "cbc";
}

// Solves the feasibility problem sum_b a_b vec_b = vec0 with cbc; returns the a_b (rounded).
static std::vector<long> solve_cbc(const std::vector<BaseVector>& bvs, const SPoly& vec0, bool msg) {
    // monomial -> constraint row
    std::vector<Exps> rows;
    std::unordered_map<Exps, int, ExpsHash> rowidx;
    for (const BaseVector& b : bvs)
        for (const auto& kv : b.vec)
            if (rowidx.emplace(kv.first, (int)rows.size()).second) rows.push_back(kv.first);

    char lp_path[] = "/tmp/schubmult_posXXXXXX.lp";
    int fd = mkstemps(lp_path, 3);
    if (fd < 0) fail("cannot create temp file");
    close(fd);
    std::string sol_path = std::string(lp_path) + ".sol";
    {
        std::ofstream lp(lp_path);
        lp << "Minimize\n obj: 0 a0\nSubject To\n";
        std::vector<std::string> row_terms(rows.size());
        for (size_t b = 0; b < bvs.size(); ++b)
            for (const auto& kv : bvs[b].vec) {
                std::string& s = row_terms[rowidx[kv.first]];
                long c = kv.second;
                if (c >= 0)
                    s += " + " + std::to_string(c) + " a" + std::to_string(b);
                else
                    s += " - " + std::to_string(-c) + " a" + std::to_string(b);
            }
        for (size_t r = 0; r < rows.size(); ++r) {
            auto it = vec0.find(rows[r]);
            if (it == vec0.end()) fail("internal: base vector monomial missing from target");
            lp << " c" << r << ":" << row_terms[r] << " = " << it->second << "\n";
        }
        lp << "Generals\n";
        for (size_t b = 0; b < bvs.size(); ++b) lp << " a" << b << "\n";
        lp << "End\n";
    }
    std::string cmd = "\"" + cbc_path() + "\" \"" + lp_path + "\" solve solu \"" + sol_path + "\"";
    if (!msg) cmd += " > /dev/null 2>&1";
    int rc = std::system(cmd.c_str());
    if (rc != 0) fail("cbc failed (set SCHUBMULT_CBC to the cbc binary)");

    std::vector<long> a(bvs.size(), 0);
    std::ifstream sol(sol_path);
    std::string line;
    if (!std::getline(sol, line)) fail("cbc produced no solution file");
    if (line.find("Infeasible") != std::string::npos || line.find("nfeasible") != std::string::npos) fail("integer program infeasible");
    while (std::getline(sol, line)) {
        std::istringstream ss(line);
        long idx;
        std::string name;
        double value;
        if (!(ss >> idx >> name >> value)) continue;
        if (name.size() > 1 && name[0] == 'a') a[std::stoul(name.substr(1))] = std::lround(value);  // round, don't truncate
    }
    std::remove(lp_path);
    std::remove(sol_path.c_str());
    return a;
}

// Y[i] = y_i, Z[i] = z_i symbols (1-indexed).
static Expr compute_positive_rep(const Expr& val, const std::vector<Expr>& Y, const std::vector<Expr>& Z, bool msg) {
    Expr val_expr = SymEngine::expand(val);
    if (SymEngine::is_a_Number(*val_expr)) return val_expr;

    // variable list: y's present sorted by index, then z's present sorted by index
    std::vector<int> yidx, zidx;
    for (const auto& s : SymEngine::free_symbols(*val_expr)) {
        const std::string& nm = SymEngine::down_cast<const SymEngine::Symbol&>(*s).get_name();
        int iy = index_of(nm, 'y'), iz = index_of(nm, 'z');
        if (iy >= 0)
            yidx.push_back(iy);
        else if (iz >= 0)
            zidx.push_back(iz);
        else
            fail("unexpected symbol " + nm);
    }
    std::sort(yidx.begin(), yidx.end());
    std::sort(zidx.begin(), zidx.end());
    int n1 = (int)yidx.size(), n3 = (int)zidx.size(), nvars = n1 + n3;
    std::map<std::string, int> varidx;
    for (int k = 0; k < n1; ++k) varidx["y_" + std::to_string(yidx[k])] = k;
    for (int i = 0; i < n3; ++i) varidx["z_" + std::to_string(zidx[i])] = n1 + i;
    int z1pos = -1;
    for (int i = 0; i < n3; ++i)
        if (zidx[i] == 1) z1pos = n1 + i;

    SPoly full = to_spoly(val_expr, varidx, nvars);
    SPoly vec0 = kill_z1(full, z1pos);

    // lookup[z-part] = monomials with that z-part whose y-part is 0/1-valued with some 1;
    // mn1L = monomials with zero y-part
    std::unordered_map<Exps, std::vector<Exps>, ExpsHash> lookup;
    std::vector<Exps> mn1L;
    for (const auto& kv : full) {
        const Exps& mm0 = kv.first;
        Exps key(mm0.begin() + n1, mm0.end());
        bool zero_one = true, has_one = false, all_zero = true;
        for (int k = 0; k < n1; ++k) {
            if (mm0[k] != 0 && mm0[k] != 1) zero_one = false;
            if (mm0[k] == 1) has_one = true;
            if (mm0[k] != 0) all_zero = false;
        }
        if (zero_one && has_one) lookup[key].push_back(mm0);
        if (all_zero) mn1L.push_back(mm0);
    }

    std::vector<BaseVector> base_vectors;
    std::map<std::vector<std::pair<int, int>>, size_t> seen;  // Python keys base vectors by the (canonicalized) product
    for (const Exps& mn1 : mn1L) {
        std::vector<std::vector<std::pair<int, int>>> comb(1);  // start with the empty product
        for (int i = n1; i < nvars; ++i) {
            if (mn1[i] == 0) continue;
            Exps mn1_2(mn1.begin() + n1, mn1.end());
            mn1_2[i - n1] = 0;
            std::vector<std::vector<std::pair<int, int>>> comb2;
            auto it = lookup.find(mn1_2);
            if (it != lookup.end())
                for (const Exps& mm0 : it->second) {
                    std::vector<std::pair<int, int>> prd;
                    for (int k = 0; k < n1; ++k)
                        if (mm0[k] == 1) prd.push_back({k, i});
                    for (const auto& a : comb) {
                        std::vector<std::pair<int, int>> f = a;
                        f.insert(f.end(), prd.begin(), prd.end());
                        comb2.push_back(std::move(f));
                    }
                }
            comb.swap(comb2);
        }
        for (auto& factors : comb) {
            std::sort(factors.begin(), factors.end());
            if (seen.count(factors)) continue;
            SPoly P;
            P[Exps(nvars, 0)] = 1;
            for (const auto& f : factors) P = mul_root(P, f.first, f.second);
            SPoly dct2 = kill_z1(P, z1pos);
            bool bad = false;
            for (const auto& kv : dct2) {
                auto it = vec0.find(kv.first);
                long have = it == vec0.end() ? 0 : it->second;
                if (std::labs(have) < std::labs(kv.second)) {
                    bad = true;
                    break;
                }
            }
            if (bad) continue;
            seen[factors] = base_vectors.size();
            base_vectors.push_back({factors, std::move(dct2)});
        }
    }
    if (base_vectors.empty()) fail("no candidate root products for " + val_expr->__str__());

    std::vector<long> a = solve_cbc(base_vectors, vec0, msg);

    SymEngine::vec_basic terms;
    for (size_t b = 0; b < base_vectors.size(); ++b) {
        if (a[b] == 0) continue;
        SymEngine::vec_basic fs;
        fs.push_back(SymEngine::integer(a[b]));
        for (const auto& f : base_vectors[b].factors) fs.push_back(SymEngine::sub(Y.at(yidx[f.first]), Z.at(zidx[f.second - n1])));
        terms.push_back(SymEngine::mul(fs));
    }
    Expr val2 = SymEngine::add(terms);
    if (!SymEngine::eq(*SymEngine::expand(SymEngine::sub(val, val2)), *SymEngine::zero)) fail("representation check failed for " + val_expr->__str__());
    return val2;
}

// schubmult_double.sv_posify (same-variable case): substitute y_i -> y_1 + r_1 + ... + r_{i-1},
// expand (the coefficient is a polynomial in the r_j with nonnegative coefficients; y_1 drops
// out), then write each monomial c * prod r_j^e_j as c * prod (y_{j+1} - y_j)^e_j.
static Expr sv_posify(const Expr& val, const std::vector<Expr>& Y) {
    Expr val_expr = SymEngine::expand(val);
    if (SymEngine::is_a_Number(*val_expr)) return val_expr;

    int maxy = 1;
    for (const auto& s : SymEngine::free_symbols(*val_expr)) {
        const std::string& nm = SymEngine::down_cast<const SymEngine::Symbol&>(*s).get_name();
        int iy = index_of(nm, 'y');
        if (iy < 0) fail("sv_posify: unexpected symbol " + nm);
        maxy = std::max(maxy, iy);
    }
    std::vector<Expr> R(maxy + 1);
    for (int j = 1; j <= maxy; ++j) R[j] = SymEngine::symbol("r_" + std::to_string(j));
    SymEngine::map_basic_basic subs;
    SymEngine::vec_basic acc{Y.at(1)};
    for (int i = 2; i <= maxy; ++i) {
        acc.push_back(R[i - 1]);
        subs[Y.at(i)] = SymEngine::add(acc);
    }
    Expr sub = SymEngine::expand(val_expr->subs(subs));
    if (SymEngine::is_a_Number(*sub)) return sub;

    std::map<std::string, int> varidx;
    for (int j = 1; j < maxy; ++j) varidx["r_" + std::to_string(j)] = j - 1;  // y_1 must be gone: to_spoly fails otherwise
    SPoly P = to_spoly(sub, varidx, std::max(maxy - 1, 1));

    // sympy's simplify pulls out the common monomial factor; do the same so r_1*r_2 + r_1**2
    // prints as (y_3 - y_2)*(y_4 - y_2) rather than as two terms.
    Exps g;
    for (const auto& kv : P) {
        if (g.empty())
            g = kv.first;
        else
            for (size_t j = 0; j < g.size(); ++j) g[j] = std::min(g[j], kv.first[j]);
    }
    auto root_pow = [&](int j, int e) { return SymEngine::pow(SymEngine::sub(Y.at(j + 1), Y.at(j)), SymEngine::integer(e)); };
    SymEngine::vec_basic terms;
    for (const auto& kv : P) {
        SymEngine::vec_basic fs;
        fs.push_back(SymEngine::integer(kv.second));
        for (int j = 1; j < maxy; ++j)
            if (kv.first[j - 1] - g[j - 1] != 0) fs.push_back(root_pow(j, kv.first[j - 1] - g[j - 1]));
        terms.push_back(SymEngine::mul(fs));
    }
    SymEngine::vec_basic fs{SymEngine::add(terms)};
    for (int j = 1; j < maxy; ++j)
        if (g[j - 1] != 0) fs.push_back(root_pow(j, g[j - 1]));
    return SymEngine::mul(fs);
}

}  // namespace positivity
