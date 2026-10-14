"""Independent checker for Lemma C (support SAT), rebuilt from scratch (no import of engine/model.py).

Setting: E8 lines in PR #1320 order, F0 deleted -> P112.  For an open case with A2 normals {a,b,c}
(capacity 2) and cap-1 normals {d,e,f}:
  class(a) = live lines perp a and not perp b, c.   O(r) = {x in {d,e,f} : r perp x}.
For r in class(a), s in class(b), a != b, r perp s, the masked entry M[r,s] = 0 receives only
the rank-1 terms x in O(r) & O(s):   sum_x u_x[r] v_x[s] = 0.
Booleans U[x,r] = (u_x[r] != 0), V[x,r] = (v_x[r] != 0).  Clause: never exactly one x with U[x,r]&V[x,s].
Term x 'serves' class a iff U[x,.] and V[x,.] each have a true entry inside class a.
Triangle clauses: for an orthogonal triangle T inside one class, all of M[T,T] = 4I is masked and receives
only P_a (rank <= 2) and the rank-1 terms, so det(4I - Q_TT) = 0 with Q = sum_x u_x v_x^T.  Expanding,
det(4I - Q) = 64 - 16 tr Q + 4 e2(Q) - det Q, and every monomial of tr, e2, det contains a full cycle
(loop, 2-cycle or 3-cycle) of Q entries; if every cycle has a zero entry the determinant is 64 != 0
(char != 2).  So some cycle has all entries nonzero, and Q[r,s] != 0 needs some x with U[x,r]&V[x,s].
Claim (UNSAT): not every class is served by >= 2 terms.
"""
import itertools, pickle, sys
from pysat.solvers import Cadical153, Glucose4


def lines():
    L = []
    for i in range(8):
        for j in range(i + 1, 8):
            for s in (2, -2):
                v = [0] * 8; v[i] = 2; v[j] = s; L.append(tuple(v))
    for t in itertools.product((1, -1), repeat=7):
        if sum(1 for q in t if q < 0) % 2 == 0:
            L.append((1,) + t)
    return L


L = lines()
ip = lambda r, s: sum(p * q for p, q in zip(L[r], L[s])) // 4
F0 = {0, 1, 26, 27, 44, 45, 54, 55}
LIVE = [r for r in range(120) if r not in F0]
assert L[0] == (2, 2, 0, 0, 0, 0, 0, 0) and L[1] == (2, -2, 0, 0, 0, 0, 0, 0) and L[27] == (0, 0, 2, -2, 0, 0, 0, 0)

OPEN = [((0, 1, 26, 28, 37, 50), (1, 1, 2, 2, 2, 1)), ((0, 1, 26, 28, 37, 54), (1, 1, 2, 2, 2, 1)),
        ((0, 1, 28, 32, 38, 47), (1, 1, 2, 2, 1, 2)), ((0, 1, 28, 38, 62, 81), (2, 1, 1, 1, 2, 2)),
        ((0, 1, 28, 40, 61, 82), (2, 1, 1, 1, 2, 2)), ((0, 2, 15, 36, 37, 50), (2, 2, 2, 1, 1, 1)),
        ((0, 2, 15, 36, 50, 51), (2, 2, 2, 1, 1, 1)), ((0, 2, 15, 36, 50, 109), (2, 2, 2, 1, 1, 1)),
        ((0, 2, 15, 36, 108, 113), (2, 2, 2, 1, 1, 1)), ((0, 2, 15, 36, 108, 114), (2, 2, 2, 1, 1, 1)),
        ((0, 2, 15, 46, 47, 52), (2, 2, 2, 1, 1, 1)), ((0, 28, 29, 38, 40, 51), (1, 1, 1, 2, 2, 2)),
        ((0, 28, 29, 40, 89, 109), (1, 2, 1, 1, 2, 2)), ((0, 28, 32, 47, 93, 104), (1, 2, 2, 2, 1, 1)),
        ((2, 3, 16, 18, 37, 50), (1, 1, 2, 2, 2, 1)), ((2, 3, 18, 22, 38, 47), (1, 1, 2, 2, 1, 2)),
        ((2, 3, 18, 22, 47, 52), (1, 1, 2, 2, 2, 1)), ((2, 3, 18, 38, 62, 97), (2, 1, 1, 1, 2, 2)),
        ((2, 3, 18, 40, 61, 98), (2, 1, 1, 1, 2, 2)), ((2, 3, 18, 50, 61, 98), (2, 1, 1, 1, 2, 2)),
        ((2, 6, 22, 29, 42, 77), (2, 2, 1, 2, 1, 1)), ((2, 6, 22, 29, 77, 110), (2, 2, 1, 2, 1, 1)),
        ((2, 6, 22, 29, 77, 118), (2, 2, 1, 2, 1, 1))]


def check(A, c, Solver):
    two = [a for a, k in zip(A, c) if k == 2]; one = [a for a, k in zip(A, c) if k == 1]
    assert sorted(c) == [1, 1, 1, 2, 2, 2]
    # shape: A2 triple (pairwise non-orthogonal), cap-1 mutually orthogonal and orthogonal to the A2
    assert all(ip(p, q) != 0 for p, q in itertools.combinations(two, 2))
    assert all(ip(p, q) == 0 for p, q in itertools.combinations(one, 2))
    assert all(ip(p, q) == 0 for p in one for q in two)
    cls = {}
    for a in two:
        for r in LIVE:
            if r not in two and ip(r, a) == 0 and any(ip(r, b) != 0 for b in two):
                assert r not in cls; cls[r] = a
    O = {r: [x for x in one if ip(r, x) == 0] for r in cls}
    nv = 0; U = {}; V = {}
    cnt = [0]
    for r in cls:
        for x in O[r]:
            nv += 1; U[x, r] = nv; nv += 1; V[x, r] = nv
    S = Solver(); ncl = 0
    cnt[0] = nv
    def new():
        cnt[0] += 1; return cnt[0]
    for r in cls:
        for s in cls:
            if cls[r] != cls[s] and ip(r, s) == 0:
                ps = []
                for x in O[r]:
                    if x in O[s]:
                        p = new(); ps.append(p)
                        S.add_clause([-p, U[x, r]]); S.add_clause([-p, V[x, s]]); S.add_clause([p, -U[x, r], -V[x, s]])
                for i in range(len(ps)):
                    S.add_clause([-ps[i]] + [ps[j] for j in range(len(ps)) if j != i]); ncl += 1
    def poss(r, s):   # literal list: Q[r,s] possibly nonzero
        out = []
        for x in O[r]:
            if x in O[s]:
                key = ('P', x, r, s)
                if key not in PV:
                    PV[key] = new()
                    S.add_clause([-PV[key], U[x, r]]); S.add_clause([-PV[key], V[x, s]]); S.add_clause([PV[key], -U[x, r], -V[x, s]])
                out.append(PV[key])
        return out
    PV = {}
    ntri = 0
    for a in two:
        La = [r for r in cls if cls[r] == a]
        for T in itertools.combinations(La, 3):
            if not all(ip(p, q) == 0 for p, q in itertools.combinations(T, 2)):
                continue
            ntri += 1
            cycles = [[(r, r)] for r in T] + [[(r, s), (s, r)] for r, s in itertools.combinations(T, 2)]
            r, s, t = T
            cycles += [[(r, s), (s, t), (t, r)], [(r, t), (t, s), (s, r)]]
            zs = []
            for cyc in cycles:
                z = new(); zs.append(z)
                for (p, q) in cyc:
                    S.add_clause([-z] + poss(p, q))   # each entry of the cycle possibly nonzero
            S.add_clause(zs)
    srv = {}
    for x in one:
        for a in two:
            z = new(); srv[x, a] = z
            W = [r for r in cls if cls[r] == a and x in O[r]]
            S.add_clause([-z] + [U[x, r] for r in W]); S.add_clause([-z] + [V[x, r] for r in W])
    for a in two:   # at least two of the three terms serve class a
        for i in range(3):
            S.add_clause([srv[one[j], a] for j in range(3) if j != i])
    res = S.solve(); S.delete()
    return res, len(cls), (ncl, ntri)


if __name__ == '__main__':
    import os
    _open = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'engine', 'k6_open.pkl')
    if os.path.exists(_open):   # the hard-coded list must be exactly the engine's 23 open cases
        assert sorted(OPEN) == sorted((tuple(A), tuple(c)) for A, c in pickle.load(open(_open, 'rb'))), 'OPEN != k6_open.pkl'
        print('OPEN list matches ../engine/k6_open.pkl (23 cases)')
    bad = 0
    for A, c in OPEN:
        r1, n, ncl = check(A, c, Cadical153)
        r2, _, _ = check(A, c, Glucose4)
        print(A, c, 'class lines', n, 'cross pairs', ncl, 'cadical', 'SAT' if r1 else 'UNSAT', 'glucose', 'SAT' if r2 else 'UNSAT', flush=True)
        bad += bool(r1) + bool(r2)
    print('LEMMA C (each class served by >=2 terms is impossible):', 'PROVED for all 23' if bad == 0 else 'FAILED')
