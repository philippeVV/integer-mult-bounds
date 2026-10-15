"""expBF exact checks (python3.12 + numpy not needed; stdlib only except reading the pinned port bundle).
(1) E8 block on a lost-free set: exact plane decoder (G/2)[S,S] (+) I_{S^c}, cost 224 = 2*112 (block 64 = 2|S|).
(2) No mod-4 bilinear lift of the plane bases with exact (mod 4) orthogonality on free pairs (a23, rW).
(3) Exact Smith caps of odd-vector pairing configurations; real a23 terms reaching ratio r_d/rho = 2, 2.5, 3.
(4) Level-1 carry of an even-coverage rectangle family exceeds half the cost (31 > 54/2).
(5) Line-model calibration: quarter decoder C=64 (depth 2), half-integral decoder C=68 (depth 1), both exact."""
import json, sys, os
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from base import *
sys.path.insert(0, RES + '/smith')
from smith import smith_vals, v2
emap = {int(k): tuple(v) for k, v in json.load(open(RES + '/lower/e8map.json')).items()}
ip = lambda a, b: sum(x * y for x, y in zip(a, b)) // 4
Rt = [emap[q] for q in P]
def check_plane_decoder(r, terms):
    E = planes(r)
    N = {}
    for d, t in terms:
        for a, x in d.items():
            for b, y in t.items(): N[(a, b)] = N.get((a, b), 0) + x * y
    for a in range(n):
        for b in range(n):
            if a == b: assert N.get((a, a), 0) == 1, (a,)
            elif not orth(E[a], E[b]): assert N.get((a, b), 0) == 0, (a, b)
    return sum(cost(r, list(d), list(t)) for d, t in terms)
# (1)
D = json.load(open(os.path.join(HERE, 'lf32_dec.json'))); S_, r, fr = D['S'], D['r'], [tuple(v) for v in D['frame']]
assert all(dot(r[a], P[a]) == 0 and dot(r[a], r[a]) == 1 for a in range(n))
assert all(ip(f, g) == 0 for f in fr for g in fr if f != g) and all(ip(f, f) == 2 for f in fr)
E = planes(r)
lost = [(a, b) for a in S_ for b in S_ if a < b and dot(P[a], P[b]) == 0 and not orth(E[a], E[b])]
assert not lost
blk = []
for f in fr:
    c = {a: ip(Rt[a], f) for a in S_ if ip(Rt[a], f)}
    blk.append(({a: Fr(v, 4) for a, v in c.items()}, {a: Fr(v) for a, v in c.items()}))
terms = blk + [({a: Fr(1)}, {a: Fr(1)}) for a in range(n) if a not in set(S_)]
C = check_plane_decoder(r, terms); Cb = sum(cost(r, list(d), list(t)) for d, t in blk)
print('(1) lost-free |S| =', len(S_), ' E8 block cost', Cb, '(= 2|S|)', ' total exact decoder cost', C)
assert Cb == 2 * len(S_) and C == 224
# (2)
def lift_consistent(r):
    E = planes(r); F = free(r); piv = {}
    Hidx = {}
    for c in range(9):
        for d_ in range(c, 9): Hidx[(c, d_)] = 2 * n * 9 + len(Hidx)
    for a in range(n):
        for b in range(a + 1, n):
            if not F[a][b]: continue
            for i, x in enumerate(E[a]):
                for j, y in enumerate(E[b]):
                    m = 0
                    for c in range(9):
                        if y >> c & 1: m ^= 1 << ((a * 2 + i) * 9 + c)
                        if x >> c & 1: m ^= 1 << ((b * 2 + j) * 9 + c)
                        for d_ in range(9):
                            if x >> c & 1 and y >> d_ & 1: m ^= 1 << Hidx[(min(c, d_), max(c, d_))]
                    rhs = (bin(x & y).count('1') // 2) & 1
                    while m:
                        h = m.bit_length() - 1
                        if h in piv: pm, pb = piv[h]; m ^= pm; rhs ^= pb
                        else: piv[h] = (m, rhs); break
                    else:
                        if rhs: return False
    return True
for nm in ('a23', 'rW'):
    ok = lift_consistent(load_r(nm)); print('(2)', nm, 'mod-4 orthogonal lift (form I+2H, z = x + 2w) exists:', ok); assert not ok
# (3)
def dotb(a, b): return bin(a & b).count('1') & 1
caps = {}
for rho in range(1, 8):
    full = list(range(1, 1 << rho)); aff = lambda v: [x for x in full if dotb(x, v)]
    cfg = {'F/F': (full, full), 'F/A': (full, aff(1)), 'A/A(a.b=1)': (aff(1), aff(1))}
    if rho >= 2: cfg['A/A(a.b=0)'] = (aff(1), aff(2))
    caps[rho] = {k: [sum(1 for t in smith_vals([[dotb(x, y) for y in Y] for x in X], 8) if t <= dd) for dd in range(4)] for k, (X, Y) in cfg.items()}
    print('(3) rho', rho, 'r_d (d=0..3):', caps[rho])
ra = load_r('a23'); Ea = planes(ra)
for R, A, rho, dd, h in [([19, 53, 99], [5, 19, 39, 44], 3, 1, 6), ([20, 24, 36, 47, 104, 110], [8, 10, 44, 60, 75], 4, 2, 10),
                          ([8, 25, 34, 42, 44, 46, 65, 66, 96], [3, 17, 38, 65, 67, 70, 84, 87, 96, 98], 6, 2, 18)]:
    M = [[dot(x, y) for b in A for y in Ea[b]] for a in R for x in Ea[a]]
    v = smith_vals(M, dd + 1); got = sum(1 for t in v if t <= dd)
    print('(3) a23 term rho', cost(ra, R, A), 'depth', dd, 'r_d', got, 'ratio', Fr(got, rho)); assert cost(ra, R, A) == rho and got == h
# (4)
cj = json.load(open(os.path.join(HERE, 'carry_1.json')))
import itertools
A_, B_, G_ = cj['A'], cj['B'], cj['G']; k = len(G_)
def inv2(G):
    k = len(G); M = [[G[i][j] % 2 for j in range(k)] + [int(i == j) for j in range(k)] for i in range(k)]
    for c in range(k):
        p = next(i for i in range(c, k) if M[i][c]); M[c], M[p] = M[p], M[c]
        for i in range(k):
            if i != c and M[i][c]: M[i] = [x ^ y for x, y in zip(M[i], M[c])]
    return [row[k:] for row in M]
Gi = inv2(G_)
A2 = [[sum(A_[a][l] * G_[l][j] for l in range(k)) % 2 for j in range(k)] for a in range(n)]
B2 = [[sum(B_[a][l] * Gi[j][l] for l in range(k)) % 2 for j in range(k)] for a in range(n)]
rects = [([A_[a][i] for a in range(n)], [B_[a][i] for a in range(n)]) for i in range(k)] + [([A2[a][j] for a in range(n)], [B2[a][j] for a in range(n)]) for j in range(k)]
rects = [(u, v) for u, v in rects if any(u) and any(v)]
F = free(ra); cnt = [[sum(u[a] * v[b] for u, v in rects) for b in range(n)] for a in range(n)]
assert all(cnt[a][b] % 2 == 0 for a in range(n) for b in range(n) if not F[a][b])
Z = [[(cnt[a][b] // 2) % 2 if not F[a][b] else 0 for b in range(n)] for a in range(n)]
rows = []
for a in range(n):
    for x in Ea[a]:
        rows.append(sum((Z[a][b] & dot(x, y)) << (2 * b + j) for b in range(n) for j, y in enumerate(Ea[b])))
val = rank2(rows); cst = sum(cost(ra, [a for a in range(n) if u[a]], [b for b in range(n) if v[b]]) for u, v in rects)
print('(4) even-coverage family:', len(rects), 'rectangles, cost', cst, ' carry rank', val, ' ratio', Fr(val, cst)); assert 2 * val > cst
# (5)
miss = sorted(set(emap) - set(P)); frame = [emap[x] for x in miss]
cq = [[ip(Rt[a], f) for f in frame] for a in range(n)]
def check_line(terms):
    N = {}
    for d, t in terms:
        for a, x in d.items():
            for b, y in t.items(): N[(a, b)] = N.get((a, b), 0) + x * y
    assert all(N.get((a, a), 0) == 1 for a in range(n))
    assert all(N.get((a, b), 0) == 0 for a in range(n) for b in range(n) if a != b and dot(P[a], P[b]) == 1)
    return sum(pair_rank([P[a] for a in d], [P[b] for b in t]) for d, t in terms), sorted({min(v2(x) for x in d.values()) + min(v2(y) for y in t.values()) for d, t in terms})
q4 = [({a: Fr(cq[a][k], 4) for a in range(n) if cq[a][k]}, {a: Fr(cq[a][k]) for a in range(n) if cq[a][k]}) for k in range(8)]
print('(5) line incumbent: cost, e_k =', check_line(q4))
Bz = [tuple(b) for b in json.load(open(os.path.join(HERE, 'lineh_best.json')))['basis']]
# integer coordinates of each root in basis Bz
def coords(v):
    M = [[Fr(Bz[j][i]) for j in range(8)] + [Fr(v[i])] for i in range(8)]
    for c in range(8):
        p = next(i for i in range(c, 8) if M[i][c] != 0); M[c], M[p] = M[p], M[c]
        M[c] = [x / M[c][c] for x in M[c]]
        for i in range(8):
            if i != c and M[i][c] != 0: M[i] = [x - M[i][c] * y for x, y in zip(M[i], M[c])]
    out = [M[i][8] for i in range(8)]; assert all(x.denominator == 1 for x in out); return [int(x) for x in out]
Cc = [coords(v) for v in Rt]
h1 = [({a: Fr(ip(Rt[a], Bz[i]), 2) for a in range(n) if ip(Rt[a], Bz[i])}, {a: Fr(Cc[a][i]) for a in range(n) if Cc[a][i]}) for i in range(8)]
print('(5) line half-integral decoder: cost, e_k =', check_line(h1))
print('PASS_EXPBF_CHECKS')
