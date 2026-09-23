#!/usr/bin/env python3
"""Design-rank checks for Section 9 (Theorem 5, Proposition 7, Table 1) of the refreshment-designs paper.

Exact rational arithmetic (fractions.Fraction); standard library only; no data, no network, no other scripts.
Run:  python3 check_design_rank.py        (exit status 1 on any failed assertion)

Model: cell means mu(e,t) = alpha(t) + g(e) + tau(s), s = t - e + 1, normalisations tau(1) = 0, g(e_1) = 0.
X_S is the incidence matrix of the fielded support S; K = ker X_S; K_tau its tenure projection.
nu(S) := dim K_tau (the dimension of the identified set of the conditioning path).

(A) Entries {1,3,4}, T = 4 (fully observed trapezoid, stride d = 1).  Co-observed tenure pairs (1,3), (2,4),
    (1,4), (1,2) connect all four tenures, yet rank = 7 < 9 parameters and nu = 2: besides the affine direction,
    delta_tau = (0,0,1,1), delta_alpha = (0,0,-1,-1), delta_g = (0,1,1) preserves every cell mean and both
    normalisations and moves the ordinary second difference tau(3) - 2 tau(2) + tau(1).  The increment graph
    has components {1,3} and {2}; Lemma 1(ii) of the companion paper gives nu = max{d, Delta_2 - w} = 2.
(B) Entries {1,5,6} and {1,5,6,9}, T = 9: nu = 1 in both.  In the four-cohort design the last cohort has no
    follow-up (w = 0 < Delta_2 - d = 3), so a last-cohort follow-up threshold is sufficient, not necessary,
    once K > 3.
(C) Every row of Table 1 recomputed exactly, including the follow-up rows.
(D) Grid: all entry sets E within {1..10} with e_1 = 1 and 3 <= K <= 5, and T = e_K, ..., e_K + 8.  On every
    design: nu = number of components of the increment graph (Lemma 1(o)); w >= Delta_2 - d implies nu = d
    (Lemma 1(i)); nu = max{d, Delta_2 - w} when K = 3 (Lemma 1(ii)); nu <= max{d, Delta_2 - w} (Lemma 1(iv)).
    The grid also counts (a) designs on which the tenure-chain condition of the previous version holds but
    nu > d, and (b) four- and five-cohort designs with w < Delta_2 - d but nu = d.
(E) Anchors in matrix form (Theorem 5(ii)): the identified set under anchors A tau = a is tau + (K_tau cap ker A).
    On the JLPS support (entries 1, 5, 13; T = 19) anchors at tau(5), tau(9), tau(13) remove the affine
    direction only (nu drops from 4 to 3); a plateau over d + 1 = 5 consecutive tenures removes all four.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import itertools, sys

fails = 0
def check(cond, msg):
    global fails
    print(("  ok    " if cond else "  FAIL  ") + msg)
    if not cond: fails += 1

# ------------------------------------------------------------------------------------------------ algebra
def rref(rows, n):
    A = [[F(x) for x in r] for r in rows]; piv = []; r = 0
    for c in range(n):
        p = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if p is None: continue
        A[r], A[p] = A[p], A[r]; pv = A[r][c]; A[r] = [x / pv for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]; A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        piv.append(c); r += 1
        if r == len(A): break
    return A[:r], piv

def rank(rows, n): return len(rref(rows, n)[1]) if rows else 0

def kernel(rows, n):
    R, piv = rref(rows, n); basis = []
    for fc in [c for c in range(n) if c not in piv]:
        v = [F(0)] * n; v[fc] = F(1)
        for i, c in enumerate(piv): v[c] = -R[i][fc]
        basis.append(v)
    return basis

def design(cells):
    E = sorted({e for e, _ in cells}); Ts = sorted({t for _, t in cells}); S = sorted({t - e + 1 for e, t in cells})
    cols = [('a', t) for t in Ts] + [('g', e) for e in E[1:]] + [('u', s) for s in S if s != 1]
    idx = {c: i for i, c in enumerate(cols)}
    rows = []
    for e, t in sorted(cells):
        r = [0] * len(cols); r[idx[('a', t)]] = 1
        if e != E[0]: r[idx[('g', e)]] = 1
        if t - e + 1 != 1: r[idx[('u', t - e + 1)]] = 1
        rows.append(r)
    return rows, cols

def trapezoid(E, T): return {(e, t) for e in E for t in range(e, T + 1)}

def nu_tau(cells, anchors=()):
    """dim of K_tau cap ker A, anchors = list of dicts {tenure: weight} (each a linear functional of tau)"""
    rows, cols = design(cells); n = len(cols)
    extra = []
    for a in anchors:
        r = [0] * n
        for s, w in a.items():
            if s != 1: r[cols.index(('u', s))] = w
        extra.append(r)
    K = kernel(rows + extra, n)
    tcols = [i for i, c in enumerate(cols) if c[0] == 'u']
    return rank([[v[i] for i in tcols] for v in K], len(tcols)) if K else 0

def increment_components(cells):
    U = sorted({t - e + 1 for e, t in cells if (e, t + 1) in cells})
    parent = {u: u for u in U}
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for (e, t) in cells:
        if (e, t + 1) not in cells: continue
        for (e2, t2) in cells:
            if t2 == t and e2 != e and (e2, t + 1) in cells:
                a, b = find(t - e + 1), find(t - e2 + 1)
                if a != b: parent[a] = b
    return len({find(u) for u in U})

def tenure_chain_condition(cells, d):
    """the condition stated in the previous version: any two observed tenures in the same residue class mod d
       are linked by a chain of co-observed pairs (two tenures observed in the same period)"""
    S = sorted({t - e + 1 for e, t in cells}); parent = {s: s for s in S}
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    Ts = {t for _, t in cells}
    for t in Ts:
        ss = [t - e + 1 for e, tt in cells if tt == t]
        for a, b in itertools.combinations(ss, 2):
            ra, rb = find(a), find(b)
            if ra != rb: parent[ra] = rb
    return all(find(a) == find(b) for a in S for b in S if (a - b) % d == 0)

def stride(E): return reduce(gcd, [e - E[0] for e in E[1:]])

# ------------------------------------------------------------------------------------------------ (A)
print("=" * 100); print("(A) Entries {1,3,4}, T = 4\n")
cells = trapezoid([1, 3, 4], 4); rows, cols = design(cells); n = len(cols); rk = rank(rows, n)
check((len(cells), n, rk, n - rk) == (7, 9, 7, 2), f"7 cells, 9 parameters, rank 7, nullity 2  (got {len(cells)}, {n}, {rk}, {n - rk})")
v = {('a', 1): 0, ('a', 2): 0, ('a', 3): -1, ('a', 4): -1, ('g', 3): 1, ('g', 4): 1, ('u', 2): 0, ('u', 3): 1, ('u', 4): 1}
resid = [sum(F(r[i]) * v[c] for i, c in enumerate(cols)) for r in rows]
check(all(x == 0 for x in resid), "delta_tau = (0,0,1,1), delta_alpha = (0,0,-1,-1), delta_g = (0,1,1) preserves every cell mean")
affine = {('a', t): -(t - 1) for t in range(1, 5)} | {('g', e): e - 1 for e in (3, 4)} | {('u', s): s - 1 for s in (2, 3, 4)}
M = [[affine[c] for c in cols], [v[c] for c in cols]]
check(rank(M, n) == 2, "it is not a multiple of the affine direction (the two span the kernel)")
d2 = v[('u', 3)] - 2 * v[('u', 2)] + 0
check(d2 == 1, "it moves tau(3) - 2 tau(2) + tau(1) by 1: ordinary curvature is not identified although d = 1")
check(tenure_chain_condition(cells, 1), "the tenure-chain condition of the previous version holds (co-observed pairs connect all tenures)")
check(increment_components(cells) == 2 and nu_tau(cells) == 2, "increment graph has 2 components = nu (Lemma 1(o)); Lemma 1(ii): max{1, 2 - 0} = 2")

# ------------------------------------------------------------------------------------------------ (B)
print(); print("=" * 100); print("(B) Last-cohort follow-up is sufficient, not necessary, when K > 3\n")
for E in ([1, 5, 6], [1, 5, 6, 9]):
    cells = trapezoid(E, 9); rows, cols = design(cells); n = len(cols); rk = rank(rows, n)
    w = 9 - E[-1]; D2 = E[1] - E[0]; d = stride(E)
    print(f"    E = {E}, T = 9: parameters {n}, rank {rk}, nullity {n - rk}, d = {d}, w = {w}, Delta_2 - d = {D2 - d}")
    check(n - rk == 1 and nu_tau(cells) == 1, f"E = {E}: nu = 1")
check(9 - 9 < (5 - 1) - 1, "for {1,5,6,9}: w = 0 < Delta_2 - d = 3, yet nu = d = 1 (earlier cohorts supply the merges)")

# ------------------------------------------------------------------------------------------------ (C)
print(); print("=" * 100); print("(C) Table 1, exact\n")
def nu_of(E, T): return nu_tau(trapezoid(E, T))
rows_tbl = []
for E, T in [([1, 2, 3, 4, 5], 5), ([1, 2, 3, 4, 5], 6), ([1, 3, 5, 7], 7), ([1, 3, 5, 7], 9), ([1, 5, 13], 13), ([1, 5, 13], 19),
             ([1, 5, 13, 21], 21), ([1, 5, 13, 21], 23), ([1, 5], 5), ([1, 5], 9), ([1, 7], 7), ([1, 7], 12)]:
    nu = nu_of(E, T); d = stride(E); rows_tbl.append((E, T, d, nu))
    print(f"    E = {E}, T = {T}: d = {d}, nu = {nu}")
    check(nu == d, f"E = {E}, T = {T}: nu = d = {d} (no follow-up requirement)")
print()
for T in range(20, 26):
    nu = nu_of([1, 5, 13, 20], T); print(f"    E = [1, 5, 13, 20], T = {T}: d = 1, w = {T - 20}, nu = {nu}")
    rows_tbl.append(([1, 5, 13, 20], T, 1, nu))
check([nu_of([1, 5, 13, 20], T) for T in range(20, 24)] == [4, 3, 2, 1], "JLPS + entry at 20: nu = 4, 3, 2, 1 for T = 20, 21, 22, 23; nu = 1 from T = 23")
check(nu_of([1, 5, 13, 22], 25) == 1 and nu_of([1, 5, 13, 22], 24) == 2, "JLPS + entry at 22: nu = 1 from T = 25 (w = 3)")
check(nu_of([1, 5, 13, 19], 30) == 2, "JLPS + entry at 19: stride 2, nu = 2 with long follow-up")
for T in range(10, 15):
    nu = nu_of([1, 7, 10], T); print(f"    E = [1, 7, 10], T = {T}: d = 3, w = {T - 10}, nu = {nu}")
check([nu_of([1, 7, 10], T) for T in range(10, 15)] == [6, 5, 4, 3, 3], "E = {1,7,10}: nu = max{3, 6 - w}")
jl = trapezoid([1, 5, 13], 19); r_, c_ = design(jl)
check((len(jl), len(c_), rank(r_, len(c_))) == (41, 39, 35), "JLPS support: 41 cells, 39 parameters, rank 35")

# identified functionals on the JLPS support

def identified(cells, lam):
    rows, cols = design(cells); n = len(cols); K = kernel(rows, n)
    return all(sum(F(lam.get(c[1], 0)) * v[i] for i, c in enumerate(cols) if c[0] == 'u') == 0 for v in K)
check(all(identified(jl, {s - 4: 1, s: -2, s + 4: 1}) for s in range(5, 16)), "JLPS: all eleven centred lag-4 second differences (s = 5..15) are identified")
check(not any(identified(jl, {s - 1: 1, s: -2, s + 1: 1}) for s in range(2, 19)), "JLPS: no ordinary second difference is identified")
check(identified(jl, {1: 1, 2: -2, 3: 1, 5: -1, 6: 2, 7: -1}), "JLPS: the cross-class contrast Delta^2 tau(2) - Delta^2 tau(6) is identified")
check(18 - nu_tau(jl) == 14, "JLPS: the identified space of tenure functionals has dimension 18 - 4 = 14")

# ------------------------------------------------------------------------------------------------ (D)
print(); print("=" * 100); print("(D) Grid of staggered trapezoids\n")
n_des = n_o = n_i = n_ii = n_iv = 0; n_chain_bad = 0; n_K4_notnec = 0; n_K3 = 0; n_iv_strict = 0
examples_chain = []; examples_K4 = []
for K in (3, 4, 5):
    for rest in itertools.combinations(range(2, 11), K - 1):
        E = [1] + list(rest); d = stride(E); D2 = E[1] - E[0]
        for T in range(E[-1], E[-1] + 9):
            cells = trapezoid(E, T); w = T - E[-1]; nu = nu_tau(cells); n_des += 1
            if nu == increment_components(cells): n_o += 1
            if w >= D2 - d:
                n_i += (nu == d)
            else:
                n_i += 1
            if K == 3:
                n_K3 += 1; n_ii += (nu == max(d, D2 - w))
            n_iv += (nu <= max(d, D2 - w)); n_iv_strict += (nu < max(d, D2 - w))
            if tenure_chain_condition(cells, d) and nu > d:
                n_chain_bad += 1
                if len(examples_chain) < 3: examples_chain.append((E, T, d, nu))
            if K >= 4 and w < D2 - d and nu == d:
                n_K4_notnec += 1
                if len(examples_K4) < 3: examples_K4.append((E, T, d, w, nu))
print(f"    designs: {n_des}  (three-cohort: {n_K3})")
check(n_o == n_des, f"Lemma 1(o): nu = #components of the increment graph on all {n_des} designs")
check(n_i == n_des, "Lemma 1(i): w >= Delta_2 - d implies nu = d on every design")
check(n_ii == n_K3, f"Lemma 1(ii): nu = max(d, Delta_2 - w) on all {n_K3} three-cohort designs")
check(n_iv == n_des, f"Lemma 1(iv): nu <= max(d, Delta_2 - w) on every design (strict on {n_iv_strict})")
print(f"    tenure-chain condition holds but nu > d: {n_chain_bad} designs, e.g. {examples_chain}")
check(n_chain_bad > 0, "the tenure-chain condition of the previous version does not imply nu = d")
print(f"    K >= 4, w < Delta_2 - d, but nu = d: {n_K4_notnec} designs, e.g. {examples_K4}")
check(n_K4_notnec > 0, "for K >= 4 the last-cohort follow-up threshold is not necessary")

# ------------------------------------------------------------------------------------------------ (E)
print(); print("=" * 100); print("(E) Anchors in matrix form on the JLPS support\n")
check(nu_tau(jl) == 4, "no anchors: nu = 4 (affine + three periodic directions)")
check(nu_tau(jl, [{5: 1}, {9: 1}, {13: 1}]) == 3, "anchors tau(5), tau(9), tau(13) (negative-control episodes): nu = 3 -- only the affine direction is removed")
check(nu_tau(jl, [{5: 1}]) == 3, "a single level anchor at tenure 5 already removes the affine direction")
check(nu_tau(jl, [{5: 1}, {13: 1, 9: -1}]) == 3, "anchors actually supplied by the panel (level at 5, increment tau(13) - tau(9)): nu = 3")
check(nu_tau(jl, [{s: 1, s + 1: -1} for s in range(14, 18)]) == 0, "plateau tau(14) = ... = tau(18) (d + 1 = 5 consecutive tenures): nu = 0")
check(nu_tau(jl, [{s: 1, s + 1: -1} for s in range(15, 18)]) == 1, "plateau over only 4 consecutive tenures: nu = 1")

print(); print("=" * 100)
print("all checks passed" if fails == 0 else f"{fails} CHECK(S) FAILED")
sys.exit(1 if fails else 0)
