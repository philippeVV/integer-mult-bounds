#!/usr/bin/env python3
"""Half-integral decoders have T = F exactly (README section 2(d)); exact check, stdlib.

Lemma.  Let W be an integral symmetric matrix with W_aa = 2 and off-diagonal entries in {0} U odd.  Then N = W/2 is a valid
masked decoder for the labels B = J + (W mod 2) (over F2), the label window is h = rank_F2 B <= rank W + 1, and every
free pair (B_ab = 0) has W_ab odd, hence N_ab != 0: the transfer graph T equals the free graph F.
Checks:
  1. E8 P120 (120 root pairs, W = Gram): N o B = diag, T = F with free degree 56, rank_Q W = 8, h = rank_F2 B = 9.
  2. Random norm-2 vector sets in Z^{9,1} (pseudo-Euclidean, metric diag(1,..,1,-1)) with pairwise inner products in
     {0} U odd: same three facts (indefinite half-integral class).
Run: python3 -I -B halfint.py   (prints PASS_HALFINT)
"""
import itertools
import random
from fractions import Fraction as Q


def e8_roots():
    R = []
    for i, j in itertools.combinations(range(8), 2):
        for a, b in itertools.product((1, -1), repeat=2):
            v = [0] * 8; v[i], v[j] = 2 * a, 2 * b; R.append(tuple(v))
    for s in itertools.product((1, -1), repeat=8):
        if s.count(-1) % 2 == 0:
            R.append(s)
    P = []
    seen = set()
    for r in R:
        if tuple(-x for x in r) not in seen:
            P.append(r); seen.add(r)
    return P                                   # coordinates scaled by 2: true inner product = dot/4


def rank_q(M):
    A = [[Q(x) for x in row] for row in M]; r = 0
    for c in range(len(A[0])):
        p = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c] / A[r][c]; A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        r += 1
    return r


def rank_f2(M):
    basis = []
    for row in M:
        x = sum((v & 1) << j for j, v in enumerate(row))
        for b in basis:
            x = min(x, x ^ b)
        if x:
            basis.append(x)
    return len(basis)


def check(W, name):
    V = len(W)
    assert all(W[a][a] == 2 for a in range(V))
    assert all(W[a][b] == 0 or W[a][b] % 2 for a in range(V) for b in range(V) if a != b)
    B = [[(1 + W[a][b]) % 2 if a != b else 1 for b in range(V)] for a in range(V)]
    N = [[Q(W[a][b], 2) for b in range(V)] for a in range(V)]
    valid = all(N[a][a] == 1 for a in range(V)) and all(N[a][b] == 0 for a in range(V) for b in range(V) if a != b and B[a][b])
    free = [[a != b and B[a][b] == 0 for b in range(V)] for a in range(V)]
    transfer = [[a != b and N[a][b] != 0 for b in range(V)] for a in range(V)]
    TeqF = free == transfer
    rW, h = rank_q(W), rank_f2(B)
    deg = [sum(r) for r in free]
    print(f'{name}: V = {V}, decoder valid {valid}, T = F {TeqF}, free degree {min(deg)}..{max(deg)} '
          f'(mean {sum(deg) / V:.1f}), rank_Q W = {rW}, h = rank_F2 B = {h} (<= rank W + 1: {h <= rW + 1})')
    return valid and TeqF and h <= rW + 1


ok = True
P = e8_roots()
W = [[sum(x * y for x, y in zip(a, b)) // 4 for b in P] for a in P]
ok &= check(W, 'E8 P120')
assert len(P) == 120 and rank_q(W) == 8

rnd = random.Random(5)
ip = lambda x, y: sum(a * b for a, b in zip(x[:-1], y[:-1])) - x[-1] * y[-1]
norm2 = [v for v in itertools.product(range(-1, 2), repeat=9) for t in (0, 1)
         for v in [v + (t,)] if ip(v, v) == 2]
for trial in range(3):
    rnd.shuffle(norm2); S = []
    for v in norm2:
        if all(ip(v, s) == 0 or ip(v, s) % 2 for s in S) and tuple(-x for x in v) not in S:
            S.append(v)
        if len(S) == 40:
            break
    ok &= check([[ip(a, b) for b in S] for a in S], f'Z^(9,1) random set {trial}')
print('PASS_HALFINT' if ok else 'FAIL_HALFINT')
