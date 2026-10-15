"""h0-generic case: some generic Y-line has H0 in its Def <=> L affine non-constant on N: L_c = w.p_c + c0 (w != 0).
Pencil blocks are then fibers of w on N (s may join one fiber or stay alone). Generic h_b lines need w ~ p_b; then the
single line F_0 = F_b = F_b* carries {h_b, h_b*, h0} (dim 2).  Enumerate every w whose fibers can carry pencil Def > 0
(i.e. some C-part's N-portion is contained in an affine hyperplane), plus all p_b, and bound
T_Y <= max over two fibers (+s) of D + D + generic."""
import numpy as np, itertools, collections
from s1lib import *
Pm = np.array([p[b] for b in N])
def normal(cs):
    A = np.hstack([np.array([p[c] for c in cs]), np.ones((len(cs), 1))])
    A = np.vstack([A, np.append(R[s], 0.0)])          # w must lie in r_s^perp
    u, sv, vt = np.linalg.svd(A); r = (sv > 1e-9).sum()
    return r - 1, vt[r:]
cands = []
for h, cp in CPART.items():
    cs = [c for c in cp if c != s]
    if len(cs) < 8: continue
    r, ns = normal(cs)
    if r <= 7:
        assert r == 7 and len(ns) == 1, (h, r)
        cands.append(ns[0])
# include p_b directions
for b in N: cands.append(np.append(p[b], 0.0))
# dedupe normals (as w part)
W = []
for v in cands:
    w = v[:-1] if len(v) == 9 else v[:8]
    w = v[:8] if len(v) == 9 else v
    w = w / np.linalg.norm(w)
    if not any(abs(abs(w @ x) - 1) < 1e-9 for x in W): W.append(w)
print('candidate directions', len(W))
best = (0, None)
for w in W:
    vals = {c: round(float(p[c] @ w), 7) for c in N}
    fib = collections.defaultdict(list)
    for c, v in vals.items(): fib[v].append(c)
    F = list(fib.values())
    gen = 2 if any(abs(abs(w @ p[b]) / np.linalg.norm(p[b]) - 1) < 1e-9 for b in N) else 1
    opts = []
    for i, f in enumerate(F): opts.append(dimD(f)); opts.append(dimD(f + [s]))
    sc = [dimD([s])]
    tot = 0
    for i, j in itertools.permutations(range(len(F)), 2):
        tot = max(tot, dimD(F[i] + [s]) + dimD(F[j]), dimD(F[i]) + dimD(F[j]))
    for i in range(len(F)): tot = max(tot, dimD(F[i]) + dimD([s]))
    if tot + gen > best[0]: best = (tot + gen, (sorted(len(f) for f in F), gen))
print('max T_Y in h0-generic case', best)
