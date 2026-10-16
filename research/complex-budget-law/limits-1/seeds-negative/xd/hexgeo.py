"""expXD: split Cayley hexagon H(q) from Q(6,q) (Van Maldeghem 2.4.13 Grassmann equations), with verification.
Points: all points of Q(6,q): X0X4+X1X5+X2X6 = X3^2.
Lines: totally singular lines of Q(6,q) whose Grassmann coords satisfy
  p12=p34, p54=p32, p20=p35, p65=p30, p01=p36, p46=p31.
Field: GF from expNC/gf.py (tables). Elements are ints."""
import sys, itertools, collections, numpy as np
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # gf.py: copy of expNC/gf.py
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
import gf
gf.IRR.setdefault((3, 3), [1, 2, 0])      # x^3 + 2x + 1
gf.IRR.setdefault((3, 4), [2, 0, 0, 1])   # x^4 + x^3 + 2 (primitive over F3)
gf.IRR.setdefault((2, 8), [1, 0, 1, 1, 1, 0, 0, 0])
GF = gf.GF

def normalize(F, v):
    v = list(v)
    for x in v:
        if x:
            ix = int(F.inv[x]); return tuple(int(F.mul[ix, y]) for y in v)
    return None

def build(q):
    p = {2: 2, 3: 3, 4: 2, 5: 5, 7: 7, 8: 2, 9: 3}[q]; e = {2: 1, 3: 1, 4: 2, 5: 1, 7: 1, 8: 3, 9: 2}[q]
    F = GF(p, e); A, M, N = F.add, F.mul, F.neg
    def Q(x): return int(F.sub[A[A[M[x[0], x[4]], M[x[1], x[5]]], M[x[2], x[6]]], M[x[3], x[3]]])
    pts = []
    for v in itertools.product(range(q), repeat=7):
        if any(v) and normalize(F, v) == v and Q(v) == 0: pts.append(v)
    pid = {v: i for i, v in enumerate(pts)}
    P = np.array(pts, dtype=np.int64)
    # bilinear form B(x,y)=x0y4+x4y0+x1y5+x5y1+x2y6+x6y2-2x3y3
    def Bv(x, Y):
        t = np.zeros(len(Y), dtype=np.int64)
        for (i, j) in [(0, 4), (4, 0), (1, 5), (5, 1), (2, 6), (6, 2)]:
            t = A[t, M[x[i], Y[:, j]]]
        two = int(A[1, 1])
        return F.sub[t, M[two, M[x[3], Y[:, 3]]]]
    def pl(x, y, i, j): return int(F.sub[M[x[i], y[j]], M[x[j], y[i]]])
    eqs = [((1, 2), (3, 4)), ((5, 4), (3, 2)), ((2, 0), (3, 5)), ((6, 5), (3, 0)), ((0, 1), (3, 6)), ((4, 6), (3, 1))]
    lines = {};
    for a in range(len(pts)):
        x = pts[a]; Bx = Bv(x, P)
        for b in np.nonzero(Bx == 0)[0]:
            b = int(b)
            if b <= a: continue
            y = pts[b]
            if any(pl(x, y, *u) != pl(x, y, *w) for u, w in eqs): continue
            plk = normalize(F, [pl(x, y, i, j) for i in range(7) for j in range(i + 1, 7)])
            if plk in lines: continue
            mem = set()
            for s in range(q):
                for t in range(q):
                    if s == 0 and t == 0: continue
                    mem.add(pid[normalize(F, [int(A[M[s, x[k]], M[t, y[k]]]) for k in range(7)])])
            assert len(mem) == q + 1
            lines[plk] = sorted(mem)
    L = list(lines.items())
    return dict(F=F, q=q, pts=pts, pid=pid, plk=[k for k, _ in L], lines=[m for _, m in L])

def incidence(H):
    NP = len(H['pts']); NL = len(H['lines'])
    adjP = [[] for _ in range(NP)]
    for l, mem in enumerate(H['lines']):
        for p in mem: adjP[p].append(l)
    H['adjP'] = adjP
    return NP, NL

def girth(adj_lists, n):
    """bipartite graph as adjacency list over n vertices; returns girth via BFS from each vertex"""
    best = 10 ** 9
    for s in range(n):
        dist = [-1] * n; par = [-1] * n; dist[s] = 0; dq = collections.deque([s])
        while dq:
            u = dq.popleft()
            if 2 * dist[u] + 1 >= best: break
            for w in adj_lists[u]:
                if dist[w] < 0: dist[w] = dist[u] + 1; par[w] = u; dq.append(w)
                elif w != par[u]: best = min(best, dist[u] + dist[w] + 1)
    return best

def graph(H, drop=()):
    NP, NL = incidence(H); dset = set(drop)
    adj = [[] for _ in range(NP + NL)]
    for l, mem in enumerate(H['lines']):
        for p in mem:
            if (p, l) in dset: continue
            adj[p].append(NP + l); adj[NP + l].append(p)
    return adj

def diameter_bip(adj, n, samples=None):
    D = 0
    for s in (range(n) if samples is None else samples):
        dist = [-1] * n; dist[s] = 0; dq = collections.deque([s])
        while dq:
            u = dq.popleft()
            for w in adj[u]:
                if dist[w] < 0: dist[w] = dist[u] + 1; dq.append(w)
        D = max(D, max(dist))
    return D

if __name__ == "__main__":
    for q in map(int, sys.argv[1:]):
        H = build(q); NP, NL = incidence(H)
        degP = collections.Counter(len(a) for a in H['adjP']); degL = collections.Counter(len(m) for m in H['lines'])
        adj = graph(H); n = NP + NL
        g = girth(adj, n); D = diameter_bip(adj, n, samples=range(0, n, max(1, n // 50)))
        exp = (q ** 6 - 1) // (q - 1)
        print(f"H({q}): points {NP} (expect {exp}) lines {NL} (expect {exp}) degP {dict(degP)} degL {dict(degL)} girth {g} diam(sampled) {D}", flush=True)
