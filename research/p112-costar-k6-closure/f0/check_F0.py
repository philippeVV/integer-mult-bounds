"""Independent checker: no nine-term masked P112 co-star decoder whose normal set contains the frame F0.
(Rebuilt from scratch; imports neither model.py, blocks.py, rsat.py nor xsat.py.)

Cases: (A) normals = F0, capacities (2,1^7) with the capacity-2 normal g0 in F0 (all 8 choices checked);
       (B) normals = F0 ∪ {x}, all capacities 1, x any of the 112 live lines (all checked).

TRIANGLE LEMMA (any field, char != 2).  Let Q1,Q3,Q5 be orthogonal 4-cliques of live lines; all lines of a
clique have the same F0-seer set Γ(Q) (F0 normals orthogonal to the line); the three cliques pairwise share
exactly the same F0 pair P; and there are 2-subsets A⊂Q1, C⊂Q3, B⊂Q5 with A∪C∪B and Ā∪C̄∪B̄ orthogonal 6-sets.
Suppose (i) Q1,Q3 are pure-like: sum_{r∈Q} v_r[P] u_r[P]^T = 4 I_2 (u_r[P] = the P-coordinates of row r of the
factor U); (ii) M[Q5,Q5] = 4I receives, besides the two P-terms, at most 3 further rank-1 terms, so
det(4I_2 - sum_{r∈Q5} v_r[P] u_r[P]^T) = 0; (iii) on the six matched part pairs, in both directions, M[r,s] = 0
receives only the P-terms: u_r[P].v_s[P] = 0.  Then contradiction.  Certificate: triangle_cert.json,
sum g_i f_i = 1 over QQ, all denominators powers of 2 (verified here by exact polynomial arithmetic).
Proof sketch: E=V1[A]^T U1[A]/4 is idempotent, the C-part of Q3 equals I-E, the B-part of Q5 equals (I-E)Π,
which forces Π = V5[P]^T U5[P]/4 = 0, so 4I - U5[P]V5[P]^T is invertible and needs 4 > 3 further terms.

Case A: Q1,Q3 not seen by g0 (seers all cap 1: U_Q square, so V^T U = 4I gives (i)); Q5 seen by g0
        (other terms: one cap-1 seer + g0 with cap 2 = 3).  Cross entries between cliques receive only the
        common seers = P.  One triangle per g0 suffices.
Case B: x has one term on S_x = {r live: r ⊥ x or r = x}.  Booleans a_r = (u_x[r] != 0), b_r = (v_x[r] != 0),
        false off S_x.  A clique is pure-like if all a_r or all b_r on it are false (then the x-term vanishes on
        M[Q,Q] and U_Q V_Q^T = 4I with square U_Q).  Q5 always has <= 2 + 1 = 3 further terms.  Clauses:
        FR: cliques with no common F0 seer are completely orthogonal and only x sees both: not(a_r & b_s).
        TRI: for each triangle and each choice of the deficient clique, not(both others pure-like and no
        a_r & b_s on the 12 directed matched part blocks).  UNSAT for every x with CaDiCaL and Glucose.
"""
import itertools, json, sys, time
from fractions import Fraction
from pysat.solvers import Cadical153, Glucose4

# ---------------- E8 lines, #1320 order ----------------
L = []
for i in range(8):
    for j in range(i + 1, 8):
        for s in (2, -2):
            v = [0] * 8; v[i] = 2; v[j] = s; L.append(tuple(v))
for t in itertools.product((1, -1), repeat=7):
    if sum(1 for q in t if q < 0) % 2 == 0:
        L.append((1,) + t)
assert len(L) == 120
def ip(r, s): return sum(p * q for p, q in zip(L[r], L[s])) // 4
def orth(r, s): return r != s and ip(r, s) == 0
F0 = (0, 1, 26, 27, 44, 45, 54, 55)
assert all(orth(a, b) for a in F0 for b in F0 if a != b)
LIVE = [r for r in range(120) if r not in F0]
SEER = {r: frozenset(f for f in F0 if ip(f, r) == 0) for r in LIVE}
assert all(len(SEER[r]) == 4 for r in LIVE)

# ---------------- cliques ----------------
blocks = {}
for r in LIVE: blocks.setdefault(SEER[r], []).append(r)
CL = []
for sv, ls in blocks.items():
    assert len(ls) == 8
    a = ls[0]; Q = [a] + [b for b in ls[1:] if orth(a, b)]
    Q2 = [b for b in ls if b not in Q]
    for X in (Q, Q2):
        assert len(X) == 4 and all(orth(p, q) for p in X for q in X if p != q)
        CL.append((sv, tuple(sorted(X))))
assert len(CL) == 28

def triangles():
    out = []
    for i, j, k in itertools.combinations(range(28), 3):
        P = CL[i][0] & CL[j][0]
        if len(P) != 2 or CL[i][0] & CL[k][0] != P or CL[j][0] & CL[k][0] != P: continue
        for A in itertools.combinations(CL[i][1], 2):
            for C in itertools.combinations(CL[j][1], 2):
                for B in itertools.combinations(CL[k][1], 2):
                    S6 = A + C + B
                    if not all(orth(p, q) for p in S6 for q in S6 if p != q): continue
                    Ab = tuple(r for r in CL[i][1] if r not in A); Cb = tuple(r for r in CL[j][1] if r not in C)
                    Bb = tuple(r for r in CL[k][1] if r not in B); T6 = Ab + Cb + Bb
                    assert all(orth(p, q) for p in T6 for q in T6 if p != q)
                    out.append((P, (i, A, Ab), (j, C, Cb), (k, B, Bb)))
    return out
TRI = triangles()

# ---------------- certificate for the abstract triangle system ----------------
def padd(p, q, c=1):
    r = dict(p)
    for e, a in q.items():
        r[e] = r.get(e, 0) + c * a
        if r[e] == 0: del r[e]
    return r
def pmul(p, q):
    r = {}
    for e1, a in p.items():
        for e2, b in q.items():
            e = tuple(x + y for x, y in zip(e1, e2)); r[e] = r.get(e, 0) + a * b
    return {e: a for e, a in r.items() if a != 0}

def check_certificate(path=__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), 'triangle_cert.json')):
    d = json.load(open(path)); names = d['names']; n = len(names); idx = {s: i for i, s in enumerate(names)}
    def var(s):
        e = [0] * n; e[idx[s]] = 1; return {tuple(e): Fraction(1)}
    const = lambda c: {tuple([0] * n): Fraction(c)} if c else {}
    ROWS = {'1': ['a1', 'a2', 'a3', 'a4'], '3': ['c1', 'c2', 'c3', 'c4'], '5': ['b1', 'b2', 'b3', 'b4']}
    PARTS = {'1': (['a1', 'a2'], ['a3', 'a4']), '3': (['c1', 'c2'], ['c3', 'c4']), '5': (['b1', 'b2'], ['b3', 'b4'])}
    assert sorted(names) == sorted(f'{t}{r}_{k}' for q in '135' for r in ROWS[q] for t in 'uv' for k in range(2))
    u = lambda r, k: var(f'u{r}_{k}'); v = lambda r, k: var(f'v{r}_{k}')
    def dot(p, q):
        s = {}
        for k in range(2): s = padd(s, pmul(p(k), q(k)))
        return s
    E = []
    for q in '13':
        for g in range(2):
            for h in range(2):
                s = const(-4 if g == h else 0)
                for r in ROWS[q]: s = padd(s, pmul(v(r, g), u(r, h)))
                E.append(s)
    N = [[None] * 2 for _ in range(2)]
    for g in range(2):
        for h in range(2):
            s = {}
            for r in ROWS['5']: s = padd(s, pmul(v(r, g), u(r, h)))
            N[g][h] = s
    E.append(padd(pmul(padd(const(4), N[0][0], -1), padd(const(4), N[1][1], -1)), pmul(N[0][1], N[1][0]), -1))
    for x, y in (('1', '3'), ('1', '5'), ('3', '5')):
        for i in range(2):
            for r in PARTS[x][i]:
                for s in PARTS[y][i]:
                    E.append(dot(lambda k: u(r, k), lambda k: v(s, k)))
                    E.append(dot(lambda k: u(s, k), lambda k: v(r, k)))
    co = d['cofactors']; assert len(co) == len(E) == 57
    tot = {}; dens = set()
    for c, f in zip(co, E):
        g = {}
        for e, q in c:
            fr = Fraction(q); g[tuple(e)] = fr; dens.add(fr.denominator)
        tot = padd(tot, pmul(g, f))
    assert tot == const(1), 'certificate does not reduce to 1'
    assert all(dd & (dd - 1) == 0 for dd in dens)
    return sum(len(c) for c in co), max(dens)

# ---------------- Case A ----------------
def case_A():
    res = {}
    for g0 in F0:
        found = None
        for t in TRI:
            P, *parts = t
            ids = [p[0] for p in parts]
            for k in range(3):
                q5 = ids[k]; others = [ids[m] for m in range(3) if m != k]
                if g0 in CL[q5][0] and all(g0 not in CL[o][0] for o in others):
                    found = (t, q5); break
            if found: break
        assert found, g0
        (P, *parts), q5 = found
        assert g0 not in P and len(CL[q5][0] - P) == 2   # Q5 extra terms: 1 (cap 1) + 2 (g0) = 3
        res[g0] = (sorted(P), [(CL[i][1], X) for i, X, _ in parts], CL[q5][1])
    return res

# ---------------- Case B ----------------
def case_B(x, Solver):
    Sx = set(r for r in LIVE if ip(x, r) == 0 or r == x)
    nv = [0]
    def new(): nv[0] += 1; return nv[0]
    a = {r: new() for r in LIVE}; b = {r: new() for r in LIVE}
    S = Solver(); add = S.add_clause
    for r in LIVE:
        if r not in Sx: add([-a[r]]); add([-b[r]])
    for i in range(28):
        for j in range(28):
            if i != j and not (CL[i][0] & CL[j][0]):
                for r in CL[i][1]:
                    for s in CL[j][1]:
                        assert orth(r, s); add([-a[r], -b[s]])
    pl = []
    for sv, Q in CL:
        na, nb, p = new(), new(), new()
        for r in Q: add([-na, -a[r]]); add([-nb, -b[r]])
        add([na] + [a[r] for r in Q]); add([nb] + [b[r] for r in Q])
        add([-p, na, nb]); add([p, -na]); add([p, -nb]); pl.append(p)
    ab = {}
    def AB(r, s):
        if (r, s) not in ab:
            t = new(); add([-t, a[r]]); add([-t, b[s]]); add([t, -a[r], -b[s]]); ab[r, s] = t
        return ab[r, s]
    for P, *parts in TRI:
        xs = []
        for (i, X, Xb), (j, Y, Yb) in itertools.combinations(parts, 2):
            for U, W in ((X, Y), (Y, X), (Xb, Yb), (Yb, Xb)):
                xs += [AB(r, s) for r in U for s in W]
        ids = [p[0] for p in parts]
        for k in range(3):
            o = [ids[m] for m in range(3) if m != k]
            add([-pl[o[0]], -pl[o[1]]] + xs)
    r = S.solve(); S.delete(); return r

if __name__ == '__main__':
    t0 = time.time()
    print('cliques 28, triangles with orthogonal 6-sets:', len(TRI))
    nt, md = check_certificate()
    print(f'triangle certificate: sum g_i f_i = 1 verified ({nt} cofactor terms, max denominator {md}, all powers of 2)')
    for g0, (P, parts, q5) in case_A().items():
        print(f'case A g0={g0}: P={P} triangle {parts} deficient clique {q5}')
    print('CASE A (normals F0, caps (2,1^7), any g0): infeasible')
    bad = []
    for x in LIVE:
        r1 = case_B(x, Cadical153); r2 = case_B(x, Glucose4)
        if r1 or r2: bad.append(x)
    print(f'case B: x over all {len(LIVE)} live lines: cadical+glucose UNSAT for {len(LIVE) - len(bad)}; SAT for {bad}')
    print('CASE B (normals F0 ∪ {x}, caps 1^9):', 'infeasible for every x' if not bad else 'OPEN')
    print('%.1fs' % (time.time() - t0))
