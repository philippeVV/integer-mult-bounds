"""Exact UPPER bound on rho'_K: every support (W', Q') of rank > c has >= c+1 rows and cols, so it is reached by removing
<= |W|-(c+1) sinks and <= |Q|-(c+1) labels.  We check every such support (up to the block stabilizer in G) for an
idempotent (rhoprime.feasible, exact over Q).  None feasible  =>  rho'_K <= c.
usage: rpub.py n orbit c"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import sys, time, itertools, core6, numpy as np
from rhoprime import feasible, rkQ
n, k, c = map(int, sys.argv[1:4])
D = core6.load(n); N = D['N']; korb = sorted(set(D['corb']))
i = D['corb'].index(korb[k]); Wm, Qm = D['cols'][i]; W = core6.bits(Wm); Q = core6.bits(Qm)
# group elements
gens = D['gens']; G = {tuple(range(N))}; fr = list(G)
while fr:
    nf = []
    for p in fr:
        for g in gens:
            q = tuple(g[p[j]] for j in range(N))
            if q not in G: G.add(q); nf.append(q)
    fr = nf
Ws, Qs = set(W), set(Q)
stab = [g for g in G if {g[w] for w in W} == Ws and {g[q] for q in Q} == Qs]
print('orbit %d beta=%d |W|=%d |Q|=%d rank=%d |stab|=%d target rho<=%d' % (k, D['beta'][i], len(W), len(Q), rkQ(W, Q), len(stab), c), flush=True)
rw, rq = len(W) - (c + 1), len(Q) - (c + 1)
t0 = time.time(); seen = set(); checked = 0; feas = []
B = {(w, q): core6.dot(w, q) for w in W for q in Q}
Bfull = np.array([[B[w, q] for q in Q] for w in W], float)
for a in range(rw + 1):
    for dw in itertools.combinations(W, a):
        kw = min(tuple(sorted(g[x] for x in dw)) for g in stab)
        if kw in seen: continue
        seen.add(kw)
        W2 = [w for w in W if w not in dw]; Wi = [W.index(w) for w in W2]
        if np.linalg.matrix_rank(Bfull[Wi]) <= c: continue
        cands = []; Mlist = []
        for b in range(rq + 1):
            for dq in itertools.combinations(range(len(Q)), b):
                keep = [j for j in range(len(Q)) if j not in dq]
                cands.append(keep)
        # batch rank by groups of equal size
        from collections import defaultdict
        bys = defaultdict(list)
        for keep in cands: bys[len(keep)].append(keep)
        for L, ks in bys.items():
            idx = np.array(ks); Ms = Bfull[Wi][:, idx].transpose(1, 0, 2)
            rr = np.linalg.matrix_rank(Ms)
            for keep, r in zip(ks, rr):
                if r <= c: continue
                # rank-maximality: if re-adding a removed sink or label keeps the rank, the larger support has the same
                # rank and (monotonicity lemma) is feasible whenever this one is; it is enumerated itself, so skip.
                Mk = Bfull[Wi][:, keep]
                if any(np.linalg.matrix_rank(np.vstack([Mk, Bfull[W.index(w)][keep]])) == r for w in dw): continue
                if any(np.linalg.matrix_rank(np.hstack([Mk, Bfull[Wi][:, [j]]])) == r for j in range(len(Q)) if j not in keep): continue
                checked += 1
                Q2 = [Q[j] for j in keep]
                if feasible(W2, Q2, exact=False) and feasible(W2, Q2): feas.append((rkQ(W2, Q2), dw, tuple(q for q in Q if q not in Q2))); print('  FEASIBLE', feas[-1], flush=True)
    print('  removals W<=%d done: W-classes %d, rank>%d checked %d (%.0f s)' % (a, len(seen), c, checked, time.time() - t0), flush=True)
print('RESULT orbit %d: %s' % (k, ('rho\' <= %d (exact: no idempotent support of rank > %d)' % (c, c)) if not feas else 'feasible supports of rank > %d exist: max %d' % (c, max(f[0] for f in feas))))
print('PASS_EXPNG_RPUB' if not feas else 'FAIL_EXPNG_RPUB')
