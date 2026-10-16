"""Cutting-plane solver for BGS dual LP (4) (unit caps unless given, d_i = 1):
   min sum cap_e w_e  s.t. y_i <= len_w(P) for every s_i-t_i path P (separated by Dijkstra), sum y_i >= 1, w >= 0.
   Options: floor=True adds w_e >= omega and, for each session, len_w(P) <= c*omega for every P in hop-shortest set H_i
   (a sufficient condition for  l = max_i dist_i <= c * min_e w_e).  mode 'maxmin': maximise omega s.t. z <= zcap."""
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
from scipy.sparse.csgraph import dijkstra
from collections import deque

def hop_paths(nv, E, s, t):
    adj = [[] for _ in range(nv)]
    for e, (u, v) in enumerate(E): adj[u].append((v, e)); adj[v].append((u, e))
    def bfs(r):
        d = [-1]*nv; d[r] = 0; dq = deque([r])
        while dq:
            x = dq.popleft()
            for y, _ in adj[x]:
                if d[y] < 0: d[y] = d[x]+1; dq.append(y)
        return d
    ds, dt = bfs(s), bfs(t); L = ds[t]; out = []
    def rec(x, p):
        if x == t: out.append(tuple(p)); return
        for y, e in adj[x]:
            if ds[y] == ds[x]+1 and dt[y] == L-ds[y]: rec(y, p+[e])
    rec(s, []); return out

def sp_path(nv, E, w, s, t):
    rows = [u for u, v in E]+[v for u, v in E]; cols = [v for u, v in E]+[u for u, v in E]
    eid = {}
    for e, (u, v) in enumerate(E): eid[(u, v)] = eid[(v, u)] = e
    G = sp.csr_matrix((np.concatenate([w, w])+1e-300, (rows, cols)), shape=(nv, nv))
    d, pred = dijkstra(G, indices=s, return_predecessors=True)
    p = []; x = t
    while x != s: p.append(eid[(pred[x], x)]); x = pred[x]
    return d[t], tuple(p)

def solve(nv, E, sess, caps=None, c=None, mode='z', zcap=None, init=None, tol=1e-8, verbose=False, method='highs'):
    m, k = len(E), len(sess); caps = np.ones(m) if caps is None else np.asarray(caps, float)
    H = [hop_paths(nv, E, s, t) for s, t in sess] if init is None else init
    cuts = [(i, p) for i in range(k) for p in H[i]]
    N = m + k + 1; om = N-1
    while True:
        R, Cc, V, b = [], [], [], []; r = 0
        for i, p in cuts:                       # y_i - len(p) <= 0
            R.append(r); Cc.append(m+i); V.append(1.)
            for e in p: R.append(r); Cc.append(e); V.append(-1.)
            b.append(0.); r += 1
        for i in range(k): R.append(r); Cc.append(m+i); V.append(-1.)
        b.append(-1.); r += 1
        if c is not None or mode == 'maxmin':
            for e in range(m): R += [r, r]; Cc += [om, e]; V += [1., -1.]; b.append(0.); r += 1
        if c is not None:
            for i in range(k):
                for p in H[i]:
                    for e in p: R.append(r); Cc.append(e); V.append(1.)
                    R.append(r); Cc.append(om); V.append(-c); b.append(0.); r += 1
        cost = np.zeros(N)
        if mode == 'maxmin':
            for e in range(m): R.append(r); Cc.append(e); V.append(caps[e])
            b.append(zcap); r += 1; cost[om] = -1
        else: cost[:m] = caps
        A = sp.csr_matrix((V, (R, Cc)), shape=(r, N))
        res = linprog(cost, A_ub=A, b_ub=np.array(b), bounds=[(0, None)]*N, method=method)
        assert res.status == 0, res.message
        w, y = res.x[:m], res.x[m:m+k]
        new = 0; dists = []
        for i, (s, t) in enumerate(sess):
            d, p = sp_path(nv, E, w, s, t); dists.append(d)
            if d < y[i]*(1-tol): cuts.append((i, p)); new += 1
        if verbose: print(f"   cp iter: obj={res.fun:.12f} new cuts={new}")
        if new == 0: return res, w, y, np.array(dists), res.x[om]
