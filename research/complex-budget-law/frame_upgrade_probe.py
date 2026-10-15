#!/usr/bin/env python3
"""Geometry probe for the 'frame-plane upgrade' experiment (self-built E8 model, not the selected labels).

Upgrade an orthogonal frame S of 8 live ports to planes E_a = <q_a, r> with one common r.
Masked-identity target (Discussion 1587): N_ab = 0 whenever E_a and E_b are not (binary) orthogonal.
The incumbent line decoder output is N = G/2. New zero constraints appear on pairs (a in S, b not in S)
with G_ab odd and r.q_b = 1. We report the size/rank of that correction block per frame.
Binary label model: q = (1, a mod 2E8) with form x0 y0 + beta, beta(x,y) = <x,y> mod 2 on E8/2E8,
so q_a.q_b = 1 + G_ab (mod 2) as in the repo's S9 = (G+1) mod 2 convention.
"""
import itertools
from fractions import Fraction as Q

roots = []
for i, j in itertools.combinations(range(8), 2):
    for si in (1, -1):
        for sj in (1, -1):
            v = [0] * 8
            v[i], v[j] = 2 * si, 2 * sj
            roots.append(tuple(v))                     # doubled coordinates
for signs in itertools.product((1, -1), repeat=8):
    if signs.count(-1) % 2 == 0:
        roots.append(tuple(signs))
assert len(roots) == 240
dot = lambda u, v: sum(a * b for a, b in zip(u, v)) // 4   # true inner product
lines = []
for v in roots:
    if tuple(-x for x in v) not in lines:
        lines.append(v)
assert len(lines) == 120
G = [[dot(u, v) for v in lines] for u in lines]

# mod-2 image: coordinates of v in a Z-basis of E8, reduced mod 2
basis = [(2, -2, 0, 0, 0, 0, 0, 0), (0, 2, -2, 0, 0, 0, 0, 0), (0, 0, 2, -2, 0, 0, 0, 0), (0, 0, 0, 2, -2, 0, 0, 0),
         (0, 0, 0, 0, 2, -2, 0, 0), (0, 0, 0, 0, 0, 2, -2, 0), (0, 0, 0, 0, 0, 2, 2, 0), (-1,) * 8]


def solve(v):
    # rational solve basis^T x = v
    n = 8
    A = [[Q(basis[j][i]) for j in range(n)] + [Q(v[i])] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        A[c] = [x / A[c][c] for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                A[r] = [x - A[r][c] * y for x, y in zip(A[r], A[c])]
    x = [A[i][n] for i in range(n)]
    assert all(t.denominator == 1 for t in x)
    return tuple(int(t) % 2 for t in x)


img = [solve(v) for v in lines]
bgram = [[dot(basis[i], basis[j]) % 2 for j in range(8)] for i in range(8)]
beta = lambda x, y: sum(x[i] * bgram[i][j] * y[j] for i in range(8) for j in range(8)) % 2
for a in range(120):
    for b in range(120):
        assert beta(img[a], img[b]) == G[a][b] % 2


def f2rank(rows):
    rows = [int(''.join(map(str, r)), 2) for r in rows if any(r)]
    rank = 0
    while rows:
        p = max(rows)
        rows = [min(r, r ^ p) for r in rows if r != p]
        rows = [r for r in rows if r]
        rank += 1
    return rank


def qrank(M):
    M = [[Q(x) for x in r] for r in M]
    rank, col = 0, 0
    rows, cols = len(M), len(M[0]) if M else 0
    for col in range(cols):
        p = next((r for r in range(rank, rows) if M[r][col] != 0), None)
        if p is None:
            continue
        M[rank], M[p] = M[p], M[rank]
        for r in range(rows):
            if r != rank and M[r][col] != 0:
                fct = M[r][col] / M[rank][col]
                M[r] = [x - fct * y for x, y in zip(M[r], M[rank])]
        rank += 1
    return rank


# frames: 8 mutually orthogonal lines
orth = [[G[a][b] == 0 for b in range(120)] for a in range(120)]
frames = set()


def extend(cur, cand):
    if len(cur) == 8:
        frames.add(tuple(sorted(cur)))
        return
    for c in cand:
        if c > cur[-1]:
            extend(cur + [c], [d for d in cand if d > c and orth[c][d]])


for a in range(120):
    extend([a], [d for d in range(a + 1, 120) if orth[a][d]])
frames = sorted(frames)
print('orthogonal frames in E8:', len(frames))
F0 = frames[0]
live = [a for a in range(120) if a not in F0]
inside = [S for S in frames if not set(S) & set(F0)]
print('frames inside P112 (disjoint from deleted frame):', len(inside))

res = []
for S in inside[:40]:
    for rp in itertools.product((0, 1), repeat=8):
        if all(beta(rp, img[a]) == 1 for a in S):
            T = [b for b in live if b not in S and beta(rp, img[b]) == 0 and any(G[a][b] % 2 for a in S)]
            block = [[G[a][b] for b in T] for a in S]
            nz = sum(1 for row in block for x in row if x % 2)
            res.append((len(T), nz, qrank(block), f2rank([[x % 2 for x in row] for row in block])))
from collections import Counter
print('per (frame, r): (|T|, new zero constraints in S x T, rank_Q, rank_F2) ->', Counter(res).most_common(8))
print('gain per upgraded frame = 8 ports * 8/9 = 7.11 deficit units per call; break-even correction charge <= 7')
