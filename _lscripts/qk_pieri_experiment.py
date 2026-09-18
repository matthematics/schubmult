"""Experiment: quantum K-theoretic top-block Pieri formula for quantum Grothendieck polynomials.

Quantum Grothendieck polynomials (Lenart--Maeno) are the FGP quantization of G_w; by
Lenart--Naito--Sagaki they represent the Schubert classes of QK(Fl_n).  We build the ring
QK(Fl_n) on the basis {G^q_w : w in S_n} from the quantum K Chevalley operators
(quantum alcove model: paths in the quantum Bruhat graph along the reduced (-omega_k)-chain,
each step weighted beta and, for a quantum edge w -> w t_ab, by q_a ... q_{b-1}), express every
class as a polynomial in the G^q_{s_k} by linear algebra, and so obtain the full product.

Sanity checks on the oracle: q = 0 recovers the K-theory products (Gx), beta = 0 recovers the
quantum cohomology products (schubmult_q), and the ring is commutative.

Then, for every u in S_n and k, expand the "top block"
    E^q_{k,k}(x; -z) G^q_u = sum_p z^{k-p} E^q_p(x_1..x_k) G^q_u,
with E^q_p = quantize(e_p(x_1..x_k)) = quantize(S_{c[p,k]}) written in the G^q basis, and compare
the coefficients with the classical rule  beta^{d-m} z^{|Fix|} (1 - beta z)^{|Left|}  (y = 0).
"""

import itertools
import sys
import time

import sympy

from schubmult import Gx, Permutation, WCGraph
from schubmult.combinatorics.permutation import uncode
from schubmult.mult.groth_double import _top_block_support, epsilon_chain
from schubmult.mult.quantum import schubmult_q
from schubmult.symbolic import sympify_sympy

n = int(sys.argv[1]) if len(sys.argv) > 1 else 3
beta = sympy.Symbol("beta")
z = sympy.Symbol("z")
q = [None] + [sympy.Symbol(f"q_{i}") for i in range(1, n)]

perms = sorted((Permutation(list(p)) for p in itertools.permutations(range(1, n + 1))), key=lambda w: (w.inv, tuple(w)))
idx = {w: i for i, w in enumerate(perms)}
N = len(perms)
ident = Permutation([])


def qk_chevalley(u, k):
    """{w: coeff of G^q_w in G^q_{s_k} * G^q_u} via quantum-Bruhat-graph paths along the (-omega_k)-chain."""
    chain = [(a, b) for a, b, _ in epsilon_chain(tuple(range(1, k + 1)), ambient_rank=n) if b <= n]
    out = {}

    def walk(start, w, size, qmon):
        if size:
            out[w] = out.get(w, 0) + beta ** (size - 1) * qmon
        for index in range(start, len(chain)):
            a, b = chain[index]
            w2 = w.swap(a - 1, b - 1)
            if w2.inv == w.inv + 1:
                walk(index + 1, w2, size + 1, qmon)
            elif w2.inv == w.inv - (2 * (b - a) - 1):
                walk(index + 1, w2, size + 1, qmon * sympy.prod([q[j] for j in range(a, b)]))

    walk(0, u, 0, sympy.Integer(1))
    return out


def chevalley_matrix(k):
    M = sympy.zeros(N, N)
    for u in perms:
        for w, c in qk_chevalley(u, k).items():
            M[idx[w], idx[u]] += c
    return M


t0 = time.time()
from sympy.polys.matrices import DomainMatrix

field = sympy.QQ.frac_field(beta, *q[1:])
Ms = [None] + [DomainMatrix.from_Matrix(chevalley_matrix(k)).convert_to(field) for k in range(1, n)]
# Artin monomials M_1^{a_1} ... M_{n-1}^{a_{n-1}}, a_k <= n - k, applied to the identity class
alphas = list(itertools.product(*[range(n - k + 1) for k in range(1, n)]))
mon_ops = {}
for alpha in alphas:
    op = DomainMatrix.eye(N, field)
    for k, a in enumerate(alpha, start=1):
        for _ in range(a):
            op = op * Ms[k]
    mon_ops[alpha] = op
id_col = idx[ident]
A = DomainMatrix.hstack(*[mon_ops[alpha][:, id_col : id_col + 1] for alpha in alphas])
A_inv = A.inv()
print(f"built Chevalley operators and inverted {N}x{N} monomial matrix in {time.time()-t0:.1f}s", flush=True)


def mult_op(x):
    """Multiplication-by-G^q_x operator as a polynomial in the Chevalley operators (DomainMatrix over Q(beta, q))."""
    rows = [[field.zero] for _ in range(N)]
    rows[idx[x]][0] = field.one
    e_x = DomainMatrix(rows, (N, 1), field)
    c = (A_inv * e_x).rep.to_list()
    T = DomainMatrix.zeros((N, N), field)
    for i, alpha in enumerate(alphas):
        if c[i][0] != field.zero:
            T = T + mon_ops[alpha] * c[i][0]
    return T.to_Matrix()


T = {x: mult_op(x) for x in perms}
print(f"multiplication operators in {time.time()-t0:.1f}s", flush=True)


def product(x, u):
    col = T[x][:, idx[u]]
    return {perms[i]: sympy.factor(col[i]) for i in range(N) if col[i] != 0}


def as_dict(elem):
    return {Permutation(list(k)): sympify_sympy(v) for k, v in elem.items() if v != 0}


# ---- sanity checks (all pairs for n <= 3, a sample otherwise) --------------------------------
import random

random.seed(0)
pairs = [(x, u) for x in perms for u in perms]
if n > 3:
    pairs = random.sample(pairs, 60)
bad_comm = bad_k = bad_qh = 0
qzero = {qi: 0 for qi in q[1:]}
b0 = {beta: 0}
sub_beta = {sympify_sympy(Gx._beta): beta}
for x, u in pairs:
    if True:
        pr = product(x, u)
        if any(sympy.simplify(T[x][idx[w], idx[u]] - T[u][idx[w], idx[x]]) != 0 for w in perms):
            bad_comm += 1
        # q = 0: classical K-theory product restricted to S_n
        ref = {w: sympify_sympy(c).subs(sub_beta) for w, c in as_dict(Gx(x) * Gx(u)).items() if w in idx}
        got = {w: c.subs(qzero) for w, c in pr.items()}
        if any(sympy.simplify(got.get(w, 0) - ref.get(w, 0)) != 0 for w in set(got) | set(ref)):
            bad_k += 1
        # beta = 0: quantum cohomology product restricted to S_n, q_j = 0 for j >= n
        refq = {}
        for w, c in schubmult_q({x: 1}, u).items():
            if w in idx:
                c = sympify_sympy(c)
                c = c.subs({sympy.Symbol(f"q_{j}"): 0 for j in range(n, 2 * n + 2)})
                if c != 0:
                    refq[w] = c
        gotq = {w: c.subs(b0) for w, c in pr.items()}
        if any(sympy.simplify(gotq.get(w, 0) - refq.get(w, 0)) != 0 for w in set(gotq) | set(refq)):
            bad_qh += 1
            if bad_qh <= 3:
                print("  QH mismatch", x, u, {tuple(w): c for w, c in gotq.items() if c != 0}, {tuple(w): c for w, c in refq.items()})
print(f"S_{n}: {len(pairs)} pairs: commutativity failures {bad_comm}, K-theory (q=0) mismatches {bad_k}, QH (beta=0) mismatches {bad_qh}", flush=True)

# ---- the top-block experiment ------------------------------------------------------------
def classical_coeff(u, w, k):
    window = [w[i] for i in range(k)]
    fixed = left = 0
    for i in range(k):
        v = u[i]
        if window[i] == v:
            fixed += 1
        elif v in window and window.index(v) < i:
            left += 1
    m = k - fixed
    d = w.inv - u.inv
    return beta ** (d - m) * z**fixed * (1 - beta * z) ** left


def top_block(u, k):
    """E^q_{k,k}(x; -z) G^q_u = sum_p z^{k-p} quantize(S_{c[p,k]}) G^q_u in the G^q basis."""
    total = {}
    for p in range(0, k + 1):
        cpk = uncode([0] * (k - p) + [1] * p)
        # S_{c[p,k]} = sum_x b_x G_x  (K-theory), quantized termwise
        for x, b in WCGraph.schub_to_groth(cpk, Gx._beta).items():
            if x not in idx:
                continue
            b = sympify_sympy(b).subs(sub_beta)
            for w, c in product(x, u).items():
                total[w] = total.get(w, 0) + z ** (k - p) * b * c
    return {w: sympy.factor(c) for w, c in total.items() if sympy.cancel(c) != 0}


def fate(u, w, k):
    window = [w[i] for i in range(k)]
    fixed = left = 0
    for i in range(k):
        v = u[i]
        if window[i] == v:
            fixed += 1
        elif v in window and window.index(v) < i:
            left += 1
    return fixed, left, k - fixed


def quantum_chains(u, k):
    """{w: [(length, q-monomial)]} over paths from u in the quantum Bruhat graph using t_ab, a <= k < b <= n,
    with b weakly decreasing along the path (the K-Pieri chain shape of elem_sym_perms_groth)."""
    out = {}

    def walk(w, last_b, length, qmon, used):
        if length:
            out.setdefault(w, []).append((length, qmon))
        for b in range(last_b, k, -1):
            for a in range(1, k + 1):
                if (a, b) in used:  # each transposition at most once (subsets of a chain); also kills 2-cycles
                    continue
                w2 = w.swap(a - 1, b - 1)
                if w2.inv == w.inv + 1:
                    walk(w2, b, length + 1, qmon, used | {(a, b)})
                elif w2.inv == w.inv - (2 * (b - a) - 1):
                    walk(w2, b, length + 1, qmon * sympy.prod([q[j] for j in range(a, b)]), used | {(a, b)})

    walk(u, n, 0, sympy.Integer(1), frozenset())
    return out


print("\n=== top block E^q_{k,k}(x;-z) G^q_u: classical-support terms vs new quantum terms ===")
nonmono = 0
for u in perms:
    for k in range(1, n):
        tb = top_block(u, k)
        supp = {w for w in _top_block_support(u, k) if w in idx}
        cls_changed = [w for w in supp if sympy.simplify(tb.get(w, 0) - classical_coeff(u, w, k)) != 0]
        missing = [w for w in supp if w not in tb]
        new_terms = {w: c for w, c in tb.items() if w not in supp}
        line = f"u={list(u)} k={k}: classical terms {'unchanged' if not cls_changed else 'CHANGED at ' + str([list(w) for w in cls_changed])}"
        if missing:
            line += f" MISSING {[list(w) for w in missing]}"
        chains = quantum_chains(u, k)
        for w, c in sorted(new_terms.items(), key=lambda t: (t[0].inv, tuple(t[0]))):
            F, L, m = fate(u, w, k)
            ratio = sympy.factor(sympy.cancel(c / (z**F * (1 - beta * z) ** L)))
            mono = ratio.is_Mul or ratio.is_Pow or ratio.is_Symbol or ratio.is_Number
            if not mono:
                nonmono += 1
            qch = [(ln, qm) for ln, qm in chains.get(w, []) if qm != 1]
            line += f"\n     G_{list(w)}: {c}   = [z^{F}(1-bz)^{L}] * {ratio}   {'MONOMIAL' if mono else 'NOT MONOMIAL'};  m={m}; quantum chains (len,q): {qch}"
        print(line, flush=True)
        if cls_changed:
            for w in cls_changed:
                print(f"     {list(w)}: got {tb.get(w, 0)}   classical {sympy.factor(classical_coeff(u, w, k))}", flush=True)
print(f"\nquantum terms whose coefficient/(fate factors) is not a single monomial: {nonmono}")

# ---- support: endpoints of quantum chains that do NOT appear (analogue of 3412 for u=1243) ----
print("\n=== quantum-chain endpoints with zero coefficient (support is smaller than chain endpoints) ===")
absent = 0
for u in perms:
    for k in range(1, n):
        tb = top_block(u, k)
        ch = quantum_chains(u, k)
        for w, lst in ch.items():
            if any(qm != 1 for _, qm in lst) and w not in tb:
                absent += 1
                print(f"u={list(u)} k={k}: w={list(w)} reachable by quantum chains {lst} but coefficient 0")
print(f"absent quantum endpoints: {absent}")


# ---- Naito--Sagaki (arXiv:2211.01578, Thm 2.13; Lenart--Maeno Conj 6.7) against the oracle ----
def ns_prec(l1, l2):
    """(a,b) < (c,d) iff b > d or (b == d and a < c)."""
    (a, b), (c, d) = l1, l2
    return b > d or (b == d and a < c)


def ns_pieri_chains(u, k):
    """All k-Pieri chains (Def. 2.9) from u inside S_n: paths in QBG with labels in L_k = {a <= k < b},
    distinct labels, b weakly decreasing, and (P2): after a repeated lower index the next label is larger."""
    labels = [(a, b) for b in range(n, k, -1) for a in range(1, k + 1)]
    out = []

    def walk(w, path, qmon):
        out.append((tuple(path), w, qmon))
        for lab in labels:
            a, b = lab
            if lab in path or (path and b > path[-1][1]):
                continue
            if len(path) >= 2 and any(path[t][0] == path[-1][0] for t in range(len(path) - 1)) and not ns_prec(path[-1], lab):
                continue  # (P2) applied to s = last position, s+1 = new label
            w2 = w.swap(a - 1, b - 1)
            if w2.inv == w.inv + 1:
                walk(w2, path + [lab], qmon)
            elif w2.inv == w.inv - (2 * (b - a) - 1):
                walk(w2, path + [lab], qmon * sympy.prod([q[j] for j in range(a, b)]))

    walk(u, [], sympy.Integer(1))
    return out


def ns_markings(path, p):
    r = len(path)
    forced = set()
    t = 0
    while t < r and all(path[i][1] == path[0][1] for i in range(t + 1)) and all(path[i][0] > path[i + 1][0] for i in range(t)):
        forced.add(t)
        t += 1
    count = 0
    for M in itertools.combinations(range(r), p):
        Ms = set(M)
        if not forced <= Ms:
            continue
        ok = True
        for s in range(r):
            if s in Ms:
                if any(path[u_][0] == path[s][0] for u_ in range(s)):
                    ok = False
                    break
            elif s < r - 1 and not ns_prec(path[s], path[s + 1]):
                ok = False
                break
        if ok:
            count += 1
    return count


print("\n=== Naito--Sagaki Pieri theorem vs oracle:  G^q_u * G^q_{c[k,p]} ===")
bad_ns = tot_ns = 0
supp_mismatch = 0
for u in perms:
    for k in range(1, n):
        chains = ns_pieri_chains(u, k)
        ns_support = set()
        for p in range(0, k + 1):
            cpk = uncode([0] * (k - p) + [1] * p)
            pred = {}
            for path, w, qmon in chains:
                mk = ns_markings(path, p)
                if mk:
                    pred[w] = pred.get(w, 0) + beta ** (len(path) - p) * mk * qmon
                    ns_support.add(w)
            got = product(cpk, u)
            tot_ns += 1
            if any(sympy.simplify(got.get(w, 0) - pred.get(w, 0)) != 0 for w in set(got) | set(pred)):
                bad_ns += 1
                if bad_ns <= 3:
                    print(f"  NS mismatch u={list(u)} k={k} p={p}: got {{{', '.join(f'{list(w)}: {c}' for w, c in got.items())}}}  pred {{{', '.join(f'{list(w)}: {sympy.factor(c)}' for w, c in pred.items())}}}")
        ours = set(top_block(u, k))
        if ours != ns_support:
            supp_mismatch += 1
            print(f"  support differs u={list(u)} k={k}: ours-NS {[list(w) for w in ours - ns_support]}, NS-ours {[list(w) for w in ns_support - ours]}")
print(f"NS theorem vs oracle: {tot_ns} products, {bad_ns} mismatches; top-block support vs NS marked-chain endpoints: {supp_mismatch} differences")
