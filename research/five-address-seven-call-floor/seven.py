#!/usr/bin/env python3
"""Seven-call layer of PR #2112's five-address four-role dyadic model.

Model: F2[Lambda], Lambda = Z A + Z P + Z R (basis), four roles, paid factor
D_N(F) = I + ([F]+1) N, N a rank-one idempotent in M4(F2), free GL4(F2) mixing
(normalised away: conjugated N stay rank-one idempotents, final mixing = I by
augmentation).  Target T = diag([A+P],[A+Q],[B+P],[B+Q]).

Step 1 (exact, Fractions): frames, lattice coordinates.
Step 2 (first order, Lambda/2Lambda): label-count profiles for a 7-call word.
Step 3: projection pi: Lambda -> Z^2 with pi(R)=0 is a group hom, hence a ring hom
  F2[Lambda] -> F2[Z^2]; every R call maps to I.  A 7-call word with m_R calls gives
  a (7-m_R)-call word over {A,P,Q,B} realising pi(T).  Enumerate all such words
  exactly (all orders, all projectors consistent with first order) in F2[Z^2].
Standard library only.  Exact arithmetic throughout.
"""
import itertools, sys, json, hashlib
from fractions import Fraction as Fr

# ---------- Step 1: exact frames and coordinates ----------
def outer(u, v):
    return tuple(tuple(Fr(a) * Fr(b) for b in v) for a in u)
def madd(*ms):
    return tuple(tuple(sum(m[i][j] for m in ms) for j in range(2)) for i in range(2))
def mscale(c, m):
    return tuple(tuple(Fr(c) * x for x in row) for row in m)
def rank2(m):
    if all(x == 0 for r in m for x in r): return 0
    return 1 if m[0][0]*m[1][1]-m[0][1]*m[1][0] == 0 else 2
h = Fr(1, 2)
FR = {
    'A': outer((1, 0), (1, 0)),
    'B': outer((0, 1), (0, 1)),
    'P': outer((h, h), (1, 1)),
    'Q': outer((h, -h), (1, -1)),
    'R': outer((1, h), (1, h)),
}
for k, m in FR.items():
    assert rank2(m) == 1, k
# coordinates in basis (A,P,R)
COORD = {'A': (1, 0, 0), 'P': (0, 1, 0), 'R': (0, 0, 1), 'B': (2, 4, -4), 'Q': (3, 3, -4)}
for k, c in COORD.items():
    assert madd(mscale(c[0], FR['A']), mscale(c[1], FR['P']), mscale(c[2], FR['R'])) == FR[k], k
# A,P,R independent over Q (as vectors in Q^4)
vec = lambda m: [m[0][0], m[0][1], m[1][0], m[1][1]]
def qrank(rows):
    rows = [list(r) for r in rows]; r = 0
    for c in range(len(rows[0])):
        piv = next((i for i in range(r, len(rows)) if rows[i][c] != 0), None)
        if piv is None: continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                f = rows[i][c] / rows[r][c]
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[r])]
        r += 1
    return r
assert qrank([vec(FR['A']), vec(FR['P']), vec(FR['R'])]) == 3
add = lambda x, y: tuple(a + b for a, b in zip(x, y))
TGT = [add(COORD['A'], COORD['P']), add(COORD['A'], COORD['Q']),
       add(COORD['B'], COORD['P']), add(COORD['B'], COORD['Q'])]
assert TGT == [(1, 1, 0), (4, 3, -4), (2, 5, -4), (5, 7, -8)]

# ---------- binary 4x4 matrices as 16-bit ints (row i = bits 4i..4i+3) ----------
n = 4
def bmul(a, b):
    out = 0
    for i in range(n):
        row = (a >> (4 * i)) & 15; acc = 0
        for k in range(n):
            if row >> k & 1: acc ^= (b >> (4 * k)) & 15
        out |= acc << (4 * i)
    return out
def brank(a):
    rows = [(a >> (4 * i)) & 15 for i in range(n)]; r = 0
    for c in range(n):
        piv = next((i for i in range(r, n) if rows[i] >> c & 1), None)
        if piv is None: continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(n):
            if i != r and rows[i] >> c & 1: rows[i] ^= rows[r]
        r += 1
    return r
def btr(a): return sum((a >> (5 * i)) & 1 for i in range(n)) & 1
def diag(bits): return sum(1 << (5 * i) for i, b in enumerate(bits) if b)
PROJ = [a for a in range(1 << 16) if brank(a) == 1 and bmul(a, a) == a]
assert len(PROJ) == 120
I4 = diag((1, 1, 1, 1)); PA = diag((1, 0, 0, 1))

# ---------- Step 2: first-order profiles ----------
def f(S): r = brank(S); return r + ((r - btr(S)) & 1)
# mod-2 directions: A=(1,0,0) P=(0,1,0) R=(0,0,1) B=0 Q=(1,1,0)
# equations: S_A+S_Q=PA, S_P+S_Q=I, S_R=0 ; S_B free.
FT = [f(S) for S in range(1 << 16)]
TR = [btr(S) for S in range(1 << 16)]
def feasible(m, S):  # can m rank-one idempotents sum to S? (necessary: rank<=m, parity)
    return m >= FT[S] and (m - TR[S]) % 2 == 0
profiles = []  # (mA,mP,mQ,mB,mR, S_Q)
for total in range(0, 8):
    for mR in range(0, total + 1, 1):
        if not feasible(mR, 0): continue
        for mA in range(total - mR + 1):
            for mP in range(total - mR - mA + 1):
                for mQ in range(total - mR - mA - mP + 1):
                    mB = total - mR - mA - mP - mQ
                    for SQ in range(1 << 16):
                        if feasible(mQ, SQ) and feasible(mA, SQ ^ PA) and feasible(mP, SQ ^ I4):
                            profiles.append((total, mA, mP, mQ, mB, mR, SQ))
by_total = {}
for p in profiles: by_total.setdefault(p[0], []).append(p)
print('first-order profiles by total:', {k: len(v) for k, v in sorted(by_total.items())})
assert min(p[0] for p in profiles if p[5] >= 2) == 6  # mR=0 words: old lattice, floor 8 (#2085)
assert min(p[0] for p in profiles if p[5] == 0) == 4  # first-order alone is weak there
# A word with mR=0 lives in the old lattice ZA+ZB+ZP (#2085: floor 8). Keep only mR>=2.
seven = [p for p in by_total[7] if p[5] >= 2]
six = [p for p in by_total[6] if p[5] >= 2]
print('7-call profiles with mR>=2:', sorted({p[1:6] for p in seven}), 'SQ choices', len(seven))

# ---------- Step 3: projected group algebra F2[Z^2], pi(R)=0 ----------
pi = lambda c: (c[0], c[1])
PI = {k: pi(COORD[k]) for k in 'APQB'}
PT = [pi(t) for t in TGT]
# group-algebra element = frozenset of exponent pairs (F2 coefficients)
ONE = frozenset([(0, 0)])
def gmul_mono(x, e): return frozenset((a + e[0], b + e[1]) for a, b in x)
def gadd(x, y): return x ^ y
ZERO = frozenset()
def mat_id(): return tuple(ONE if i == j else ZERO for i in range(n) for j in range(n))
def call(M, F, N):
    """M * (I + ([F]+1) N) = M + (M N)([F]+1)."""
    e = PI[F]
    cols = [(j) for j in range(n)]
    MN = []
    for i in range(n):
        for j in range(n):
            acc = ZERO
            for k in range(n):
                if (N >> (4 * k + j)) & 1: acc = acc ^ M[i * n + k]
            MN.append(acc)
    return tuple(M[t] ^ MN[t] ^ gmul_mono(MN[t], e) for t in range(n * n))
TARGET = tuple(frozenset([PT[i]]) if i == j else ZERO for i in range(n) for j in range(n))

def tuples_summing(m, S):
    """All ordered m-tuples of projectors with XOR S (as sorted multisets -> orders later)."""
    out = set()
    def rec(k, start, acc, cur):
        if k == 0:
            if acc == S: out.add(tuple(cur))
            return
        for idx in range(start, 120):
            rec(k - 1, idx, acc ^ PROJ[idx], cur + [idx])
    rec(m, 0, 0, [])
    return out  # multisets (nondecreasing index tuples)

def words_for(mA, mP, mQ, mB, SQ):
    """Yield all distinct words (sequences of (label, projector)) with this first-order data."""
    sets = {'A': tuples_summing(mA, SQ ^ PA), 'P': tuples_summing(mP, SQ ^ I4),
            'Q': tuples_summing(mQ, SQ)}
    Bms = list(itertools.combinations_with_replacement(range(120), mB))
    for a in sets['A']:
        for p in sets['P']:
            for q in sets['Q']:
                for b in Bms:
                    items = [('A', x) for x in a] + [('P', x) for x in p] + [('Q', x) for x in q] + [('B', x) for x in b]
                    for perm in set(itertools.permutations(items)):
                        yield perm

def search(profs, label):
    checked = 0; hits = []
    seen = set()
    for (_, mA, mP, mQ, mB, mR, SQ) in profs:
        key = (mA, mP, mQ, mB, SQ)
        if key in seen: continue
        seen.add(key)
        # DFS with prefix sharing over distinct words
        words = sorted(set(words_for(mA, mP, mQ, mB, SQ)))
        cache = {(): mat_id()}
        for w in words:
            # compute product incrementally using longest cached prefix
            k = len(w)
            while w[:k] not in cache: k -= 1
            M = cache[w[:k]]
            for t in range(k, len(w)):
                M = call(M, w[t][0], PROJ[w[t][1]])
                if t < len(w) - 1: cache[w[:t + 1]] = M
            checked += 1
            if M == TARGET: hits.append(w)
            if len(cache) > 200000: cache = {(): mat_id()}
    print(f'{label}: residual profiles {len(seen)}, words checked {checked}, hits {len(hits)}')
    return checked, hits

if __name__ == '__main__':
    # positive control: direct 8-call word in projected model (roles i: two diagonal calls)
    E = [diag(tuple(int(i == j) for j in range(4))) for i in range(4)]
    direct = [('A', E[0]), ('P', E[0]), ('A', E[1]), ('Q', E[1]), ('B', E[2]), ('P', E[2]), ('B', E[3]), ('Q', E[3])]
    M = mat_id()
    for F, N in direct: M = call(M, F, N)
    assert M == TARGET, 'direct control failed'
    # negative control: drop one call
    M = mat_id()
    for F, N in direct[1:]: M = call(M, F, N)
    assert M != TARGET
    print('controls ok')
    c6, h6 = search(six, '6-call (mR>=2) projected residual')
    c7, h7 = search(seven, '7-call (mR>=2) projected residual')
    print(json.dumps({'six_checked': c6, 'six_hits': len(h6), 'seven_checked': c7, 'seven_hits': len(h7)}))
    if h7:
        print('HITS (first 5):', h7[:5])
