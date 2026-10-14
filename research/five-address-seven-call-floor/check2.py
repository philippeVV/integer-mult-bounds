#!/usr/bin/env python3
"""Independent second check of the seven-call exclusion (different enumeration,
different product algorithm).  Matrices are lists of lists of 0/1; projectors are
enumerated as outer products u v^T with v.u = 1 (not by idempotence test).

Product formula: prod_j (I + x_j N_j), x_j = [F_j] + 1 in char 2, expands as
  sum_{S increasing} (prod_{j in S} N_j) * prod_{j in S} ([F_j]+1)
  = sum_S N_S * sum_{T subset S} [sum_{j in T} F_j].
"""
import itertools, hashlib, json
from collections import Counter

V = [v for v in itertools.product((0, 1), repeat=4) if any(v)]
PRJ = set()
for u in V:
    for v in V:
        if sum(a * b for a, b in zip(u, v)) % 2 == 1:
            PRJ.add(tuple(tuple(a * b for b in v) for a in u))
PRJ = sorted(PRJ)
assert len(PRJ) == 120
Z4 = tuple(tuple(0 for _ in range(4)) for _ in range(4))
I4 = tuple(tuple(int(i == j) for j in range(4)) for i in range(4))
def mx(a, b): return tuple(tuple(x ^ y for x, y in zip(r, s)) for r, s in zip(a, b))
def mm(a, b): return tuple(tuple(sum(a[i][k] & b[k][j] for k in range(4)) & 1 for j in range(4)) for i in range(4))
def dg(*b): return tuple(tuple(b[i] if i == j else 0 for j in range(4)) for i in range(4))
PA = dg(1, 0, 0, 1)
PS = set(PRJ)

# projected exponents in Z^2 (pi(R)=0): A=(1,0) P=(0,1) B=(2,4) Q=(3,3)
E = {'A': (1, 0), 'P': (0, 1), 'B': (2, 4), 'Q': (3, 3)}
TGT = [(1, 1), (4, 3), (2, 5), (5, 7)]

def product_coeffs(word):
    """Return dict (i,l) -> set of exponents with coefficient 1."""
    k = len(word)
    acc = Counter()
    for r in range(k + 1):
        for Sset in itertools.combinations(range(k), r):
            NS = I4
            for j in Sset: NS = mm(NS, word[j][1])
            nz = [(i, l) for i in range(4) for l in range(4) if NS[i][l]]
            if not nz: continue
            for rt in range(len(Sset) + 1):
                for T in itertools.combinations(Sset, rt):
                    e = (sum(E[word[j][0]][0] for j in T), sum(E[word[j][0]][1] for j in T))
                    for il in nz: acc[(il, e)] ^= 1
    out = {}
    for (il, e), c in acc.items():
        if c: out.setdefault(il, set()).add(e)
    return out
TARGET = {((i, i)): {TGT[i]} for i in range(4)}

def is_target(word): return product_coeffs(word) == TARGET

# controls
D = [dg(*[int(i == j) for j in range(4)]) for i in range(4)]
direct = [('A', D[0]), ('P', D[0]), ('A', D[1]), ('Q', D[1]), ('B', D[2]), ('P', D[2]), ('B', D[3]), ('Q', D[3])]
assert is_target(direct) and not is_target(direct[:-1])

# old-lattice (no R) first-order floor, basis (A,B,P): Q=(1,1,1); masks
def rank(m):
    rows = [list(r) for r in m]; r = 0
    for c in range(4):
        p = next((i for i in range(r, 4) if rows[i][c]), None)
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(4):
            if i != r and rows[i][c]: rows[i] = [x ^ y for x, y in zip(rows[i], rows[r])]
        r += 1
    return r
def fl(S): r = rank(S); t = sum(S[i][i] for i in range(4)) & 1; return r + ((r - t) & 1)
ALLM = [tuple(tuple((b >> (4 * i + j)) & 1 for j in range(4)) for i in range(4)) for b in range(1 << 16)]
old = min(fl(SQ) + fl(mx(SQ, dg(1, 0, 0, 1))) + fl(mx(SQ, dg(0, 1, 1, 0))) + fl(mx(SQ, I4)) for SQ in ALLM)
print('old lattice {A,B,P,Q} first-order floor:', old)
assert old == 8

# enumerate 5-call non-R residual words independently:
# constraints (from Lambda/2Lambda): xor(A projs)+xor(Q)=PA, xor(P)+xor(Q)=I, B free; count profiles
def xor_all(ms):
    s = Z4
    for m in ms: s = mx(s, m)
    return s
def gen(counts):
    """counts: dict label->m. Enumerate all label sequences and projector tuples (ordered)
    with Q free (all but nothing determined), A,P last element determined."""
    labels = [l for l in 'APQB' for _ in range(counts[l])]
    seqs = set(itertools.permutations(labels))
    mA, mP, mQ, mB = counts['A'], counts['P'], counts['Q'], counts['B']
    for qs in itertools.product(PRJ, repeat=mQ):
        SQ = xor_all(qs)
        def det_list(m, S):
            if m == 0:
                return [()] if S == Z4 else []
            res = []
            for pre in itertools.product(PRJ, repeat=m - 1):
                last = mx(S, xor_all(pre))
                if last in PS: res.append(pre + (last,))
            return res
        As = det_list(mA, mx(SQ, PA)); Ps = det_list(mP, mx(SQ, I4))
        if not As or not Ps: continue
        for a in As:
            for p in Ps:
                for b in itertools.product(PRJ, repeat=mB):
                    for seq in seqs:
                        it = {'A': iter(a), 'P': iter(p), 'Q': iter(qs), 'B': iter(b)}
                        yield tuple((l, next(it[l])) for l in seq)

# independent profile derivation in Lambda/2Lambda (basis A,P,R): A=(1,0,0) P=(0,1,0) R=(0,0,1) B=0 Q=(1,1,0)
FL = [fl(m) for m in ALLM]; TRc = [sum(m[i][i] for i in range(4)) & 1 for m in ALLM]
idx = {m: i for i, m in enumerate(ALLM)}
def ok(m, i): return m >= FL[i] and (m - TRc[i]) % 2 == 0
iPA = idx[PA]; iI = idx[I4]
prof = set()
for total in (6, 7):
    for mR in range(total + 1):
        if not ok(mR, 0) or mR == 0: continue   # mR=0: old lattice, floor 8 above
        for mA in range(total - mR + 1):
            for mP in range(total - mR - mA + 1):
                for mQ in range(total - mR - mA - mP + 1):
                    mB = total - mR - mA - mP - mQ
                    if any(ok(mQ, q) and ok(mA, q ^ iPA) and ok(mP, q ^ iI) for q in range(1 << 16)):
                        prof.add((total, mA, mP, mQ, mB, mR))
print('profiles (total,mA,mP,mQ,mB,mR):', sorted(prof))
assert sorted(prof) == [(6, 0, 2, 2, 0, 2), (7, 0, 2, 2, 1, 2), (7, 1, 1, 3, 0, 2), (7, 1, 3, 1, 0, 2)]
results = {}
allwords = set()
for name, counts, total in [('6:(Q2,P2)', dict(A=0, P=2, Q=2, B=0), 6),
                            ('7:(Q2,P2,B1)', dict(A=0, P=2, Q=2, B=1), 7),
                            ('7:(A1,P1,Q3)', dict(A=1, P=1, Q=3, B=0), 7),
                            ('7:(A1,P3,Q1)', dict(A=1, P=3, Q=1, B=0), 7)]:
    ws = set(gen(counts))
    hits = [w for w in ws if is_target(w)]
    results[name] = (len(ws), len(hits))
    if total == 7: allwords |= ws
    print(name, 'distinct words', len(ws), 'hits', len(hits), flush=True)
print(json.dumps({'results': results, 'seven_total_distinct': len(allwords)}))
