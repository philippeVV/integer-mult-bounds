"""expXA bw16.py: Barnes-Wall BW16 (k=4) under the two label choices.

Construction D: BW16 = RM(1,4) + 2 RM(3,4) + 4 Z^16; 4320 minimal vectors of norm 8, halved Gram -> norm 4, ips {0,+-1,+-2}.
(a) parity labels B = (W+J) mod 2 (expM3/expL3): h = 9. G/4 is INVALID: ip +-2 pairs are required (B=1) but N != 0.
    We count them and also report what G/4 'would' cost (to show validity, not cost, is the obstruction).
(b) zero labels B0 = I + [ip = 0] (the least restrictive label compatible with G/4): h0 = rank_F2 B0,
    cost of the 16 coordinate terms of G/4, gamma0; greedy lower bound on the pairwise-non-orthogonal clique
    (any label compatible with G/4 has h >= that clique).
"""
import itertools, random
import numpy as np
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from common import rank_f2_rows, block_rank

n = 16
pts = list(itertools.product((0, 1), repeat=4))
rm1 = set()
for a0 in (0, 1):
    for a in itertools.product((0, 1), repeat=4):
        rm1.add(tuple((a0 + sum(ai * pi for ai, pi in zip(a, p))) % 2 for p in pts))
w8 = [c for c in rm1 if sum(c) == 8]
assert len(w8) == 30
vecs = []
for i, j in itertools.combinations(range(n), 2):
    for si in (2, -2):
        for sj in (2, -2):
            v = [0] * n; v[i] = si; v[j] = sj; vecs.append(v)
for c in w8:
    sup = [i for i in range(n) if c[i]]
    for k in range(0, 9, 2):
        for m in itertools.combinations(sup, k):
            v = [0] * n
            for i in sup: v[i] = 1
            for i in m: v[i] = -1
            vecs.append(v)
assert len(vecs) == 4320
X = np.array(vecs)
# one representative per line
seen, lines = set(), []
for v in vecs:
    t = tuple(v); nt = tuple(-x for x in v)
    if nt in seen: continue
    seen.add(t); lines.append(v)
L = np.array(lines)
V = len(L)
W = (L @ L.T) // 2
vals, cnts = np.unique(W[~np.eye(V, dtype=bool)], return_counts=True)
print("BW16 lines V =", V, " halved Gram off-diag values:", dict(zip(vals.tolist(), (cnts // 1).tolist())))

def bits_from(M):
    out = []
    for a in range(V):
        row = np.nonzero(M[a])[0]
        x = 0
        for b in row: x |= 1 << int(b)
        out.append(x)
    return out

Bpar = ((W + 1) % 2).astype(np.int8)
Bp = bits_from(Bpar)
hpar = rank_f2_rows(Bp)
bad = int(((np.abs(W) == 2) & (Bpar == 1)).sum()) // 2
print(f"(a) parity labels: h = {hpar}; pairs with ip=+-2 (required, B=1) where G/4 = +-1/2 != 0: {bad} -> G/4 INVALID")
R = [[int(x) for x in np.nonzero(L[:, j])[0]] for j in range(n)]
cp = sum(block_rank(Bp, r, r) for r in R)
print(f"    (what the 16 coordinate terms of G/4 would cost under parity labels: {cp}, 'gamma' {cp/V:.4f} -- invalid)")

B0 = ((W == 0)).astype(np.int8)  # diagonal W=4 -> set below
np.fill_diagonal(B0, 1)
Bz = bits_from(B0)
h0 = rank_f2_rows(Bz)
c0 = [block_rank(Bz, r, r) for r in R]
print(f"(b) zero labels B0 = I + [ip=0]: h0 = {h0}; G/4 valid; coordinate-term costs {sorted(set(c0))}, C = {sum(c0)}, gamma = {sum(c0)/V:.4f}")
# greedy clique in non-orthogonal graph (lower bound for every label compatible with G/4)
nonorth = (W != 0)
best = 0
rng = random.Random(1)
for trial in range(200):
    order = list(range(V)); rng.shuffle(order)
    cl = []
    cand = np.ones(V, dtype=bool)
    for a in order:
        if cand[a]:
            cl.append(a); cand &= nonorth[a]
    best = max(best, len(cl))
print(f"    greedy pairwise non-orthogonal clique: {best}  => every label compatible with G/4 has h >= {best}")
print(f"    super-frame factor V/(h*rank): parity (invalid) {V/(hpar*16):.2f}; zero labels {V/(h0*16):.3f}")

assert V == 2160 and hpar == 9 and bad == 302400 and cp == 144 and h0 == 191 and sum(c0) == 2608
print('PASS_EXPXA_BW16')
