"""expXA fixedh.py: can the Lorentzian G/2 mechanism compress at FIXED window h (multiplicity instead of growth)?

Candidates: norm-2 vectors of I_{n,1} (form diag(-1,1^n)), 0 <= x0 <= X0, |x_i| <= 2, one per line.
Valid set: every pair has W_ab = 0 or W_ab odd ("even => zero"); then N = W/2 is an exact masked decoder for the
labels B = (W+J) mod 2.  Window cap: reductions mod 2 must stay inside a span of dim <= hmax-1 (so h <= hmax).
Greedy random restarts maximise V; cost = sum over the n+1 coordinate terms of rank_F2 B[R_c,R_c] (exact).
usage: fixedh.py n hmax restarts [X0]
"""
import sys, itertools, random
import numpy as np
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from common import rank_f2_rows, block_rank

n, hmax, restarts = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
X0 = int(sys.argv[4]) if len(sys.argv) > 4 else 3
cands = []
for x0 in range(0, X0 + 1):
    target = 2 + x0 * x0
    # vectors x in {-2..2}^n with sum x_i^2 = target
    def rec(i, rem, cur):
        if i == n:
            if rem == 0:
                cands.append([x0] + cur)
            return
        for v in (0, 1, -1, 2, -2):
            if v * v <= rem:
                rec(i + 1, rem - v * v, cur + [v])
    rec(0, target, [])
# one representative per line: first nonzero entry positive
C = []
for v in cands:
    nz = next(c for c in v if c != 0)
    if nz > 0:
        C.append(v)
C = np.array(C, dtype=np.int64)
eps = np.array([-1] + [1] * n)
CE = C * eps
print(f"n={n} candidates={len(C)}", flush=True)
red = [int(''.join(str(int(b)) for b in (v % 2)), 2) for v in C]


def span_add(piv, v):
    w = v
    while w:
        lb = w.bit_length() - 1
        if lb in piv:
            w ^= piv[lb]
        else:
            return w, lb
    return 0, None


best = (0, None)
rng = random.Random(7)
for t in range(restarts):
    order = list(range(len(C))); rng.shuffle(order)
    chosen, piv = [], {}
    M = np.zeros((0, n + 1), dtype=np.int64)
    for a in order:
        w, lb = span_add(piv, red[a])
        if w and len(piv) >= hmax - 1:
            continue
        if len(chosen):
            g = M @ CE[a]
            if not np.all((g == 0) | (g % 2 != 0)):
                continue
        chosen.append(a); M = np.vstack([M, C[a]])
        if w:
            piv[lb] = w
    if len(chosen) > best[0]:
        best = (len(chosen), chosen)
V, S = best
S = sorted(S)
W = (C[S] * eps) @ C[S].T
assert all(W[i, i] == 2 for i in range(V))
assert np.all((W == 0) | (W % 2 != 0) | np.eye(V, dtype=bool))
Bb = []
for a in range(V):
    x = 0
    for b in np.nonzero((W[a] + 1) % 2)[0]:
        x |= 1 << int(b)
    Bb.append(x)
h = rank_f2_rows(Bb)
X = C[S]
cost = 0
for c in range(n + 1):
    R = [int(i) for i in np.nonzero(X[:, c])[0]]
    if R:
        cost += block_rank(Bb, R, R)
rk = np.linalg.matrix_rank(W.astype(float))
print(f"best V={V} h={h} rank(W)={rk} coordinate-term cost C={cost} gamma={cost/V:.4f} merit gamma*(h-1)={cost/V*(h-1):.3f}")
