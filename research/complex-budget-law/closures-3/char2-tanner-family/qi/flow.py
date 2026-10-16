"""Undirected multiflow F_* of the Singer PG(2,q)-M Tanner gadget (stronger-rank-screens §2 / expNJ C1).
For an information set R (|R| = S, the session arcs) each session k = R-arc (p,l) is split into two leaves;
demand k must route p -> l through the E-S internal edges (capacity 1).  F_* = W  =>  s >= Wm for every
rational frame assignment, every m (C1).
usage: flow.py e j m0 nseeds [lp 0/1]"""
import sys, random, collections, time, numpy as np
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # qb.py (expQB), gf.py (expNC)
from qb import model
e, j, m0, nseeds = map(int, sys.argv[1:5]); dolp = int(sys.argv[5]) if len(sys.argv) > 5 else 1
t0 = time.time()
Mo = model(e, j, m0, verbose=False); F = Mo['F']; N = Mo['N']; k = Mo['k']
edges = sorted(Mo['gE']); E = len(edges)
M0 = np.zeros((E, 2 * N * k), dtype=np.int16)
for t, (i, d) in enumerate(edges):
    jj = (i + d) % N
    M0[t, i * k:(i + 1) * k] = Mo['gE'][(i, d)]; M0[t, N * k + jj * k:N * k + (jj + 1) * k] = Mo['hE'][(i, d)]
K = F.nullspace(M0); S = K.shape[0]
Wb = np.zeros((S, E), dtype=np.int16)
for t, (i, d) in enumerate(edges):
    acc = np.zeros(S, dtype=np.int16)
    for r, c in enumerate(Mo['gE'][(i, d)]): acc = F.add[acc, F.mul[K[:, i * k + r], c]]
    Wb[:, t] = acc
assert F.rank(Wb) == S
print(f"q={Mo['q']} j={j} m0={m0} N={N} E={E} S={S}  6S-E={6*S-E}  distance cap (E-S)/5={(E-S)/5:.1f}  ({time.time()-t0:.1f}s)", flush=True)
ends = [(i, N + (i + d) % N) for (i, d) in edges]

def greedy(R, rng, tries=200):
    """integral edge-disjoint routing of the R-demands (random-order shortest paths); returns #routed best."""
    Rs = set(R); best = 0
    for _ in range(tries):
        used = set(); order = list(R); rng.shuffle(order); got = 0
        adj = collections.defaultdict(list)
        for t in range(E):
            if t in Rs: continue
            a, b = ends[t]; adj[a].append((b, t)); adj[b].append((a, t))
        for r in order:
            s_, g_ = ends[r]; prev = {s_: None}; dq = collections.deque([s_])
            while dq and g_ not in prev:
                x = dq.popleft()
                for y, t in adj[x]:
                    if t in used or y in prev: continue
                    prev[y] = (x, t); dq.append(y)
            if g_ in prev:
                got += 1; y = g_
                while prev[y] is not None: x, t = prev[y]; used.add(t); y = x
        best = max(best, got)
        if best == len(R): break
    return best

def pathfinder(R, rng, iters=300):
    """negotiated-congestion integral routing; returns (#demands, paths) when all edge-disjoint, else None."""
    import heapq
    Rs = set(R); adj = collections.defaultdict(list)
    for t in range(E):
        if t in Rs: continue
        a, b = ends[t]; adj[a].append((b, t)); adj[b].append((a, t))
    hist = collections.defaultdict(float); paths = {}
    for it in range(iters):
        use = collections.Counter()
        for r in paths.values(): use.update(r)
        order = list(R); rng.shuffle(order)
        for r in order:
            if r in paths: use.subtract(paths[r])
            s_, g_ = ends[r]; dist = {s_: 0.0}; prev = {}; pq = [(0.0, s_)]
            while pq:
                d0, x = heapq.heappop(pq)
                if x == g_: break
                if d0 > dist[x]: continue
                for y, t in adj[x]:
                    c = (1 + hist[t]) * (1 + 2.0 * use[t] * (it + 1))
                    if d0 + c < dist.get(y, 1e18): dist[y] = d0 + c; prev[y] = (x, t); heapq.heappush(pq, (d0 + c, y))
            p = []; y = g_
            while y != s_: x, t = prev[y]; p.append(t); y = x
            paths[r] = p; use.update(p)
        over = [t for t, c in use.items() if c > 1]
        if not over:
            # verify: edge-disjoint, each path connects its pair through internal edges
            allp = [t for p in paths.values() for t in p]; assert len(allp) == len(set(allp)) and not (set(allp) & Rs)
            for r, p in paths.items():
                cur = {ends[r][0]}; deg = collections.Counter()
                for t in p: deg[ends[t][0]] += 1; deg[ends[t][1]] += 1
                odd = {v for v, c in deg.items() if c % 2}; assert odd == {ends[r][0], ends[r][1]}
            return it, paths
        for t in over: hist[t] += 1.0
    return None

def lp(R):
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix
    Rs = set(R); I = [t for t in range(E) if t not in Rs]; nI = len(I); V = 2 * N; nk = len(R)
    # vars per commodity: 2*nI arc flows + 1 throughput
    nv = nk * (2 * nI + 1); rows = []; cols = []; vals = []; beq = []; r0 = 0
    for c, rr in enumerate(R):
        base = c * (2 * nI + 1); s_, g_ = ends[rr]
        for u, t in enumerate(I):
            a, b = ends[t]
            for (x, y, off) in ((a, b, 0), (b, a, 1)):
                v = base + 2 * u + off
                rows += [r0 + x, r0 + y]; cols += [v, v]; vals += [-1.0, 1.0]   # out of x, into y
        thr = base + 2 * nI
        rows += [r0 + s_, r0 + g_]; cols += [thr, thr]; vals += [1.0, -1.0]
        r0 += V
    Aeq = coo_matrix((vals, (rows, cols)), shape=(r0, nv)).tocsr()
    rows = []; cols = []
    for c in range(nk):
        base = c * (2 * nI + 1)
        for u in range(nI): rows += [u, u]; cols += [base + 2 * u, base + 2 * u + 1]
    Aub = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(nI, nv)).tocsr()
    cost = np.zeros(nv); bounds = [(0, None)] * nv
    for c in range(nk): cost[c * (2 * nI + 1) + 2 * nI] = -1; bounds[c * (2 * nI + 1) + 2 * nI] = (0, 1)
    res = linprog(cost, A_ub=Aub, b_ub=np.ones(nI), A_eq=Aeq, b_eq=np.zeros(r0), bounds=bounds, method='highs')
    return -res.fun

cert = 0
for seed in range(nseeds):
    rng = random.Random(seed); perm = list(range(E)); rng.shuffle(perm)
    _, piv = F.rref(Wb[:, perm]); R = sorted(perm[c] for c in piv)
    g = greedy(R, rng, 20)
    out = f"seed {seed}: |R|={len(R)} greedy integral routed={g}"
    if g < len(R):
        pf = pathfinder(R, rng)
        if pf: g = len(R); out += f"  pathfinder: ALL {len(R)} edge-disjoint (iter {pf[0]}, total length {sum(map(len, pf[1].values()))}) -> F*=W certified"; cert += 1
    if dolp and g < len(R): out += f"  LP F*={lp(R):.4f}"
    print(out + f"  ({time.time()-t0:.1f}s)", flush=True)
print('PASS_EXPQI_FLOW' if cert == nseeds else 'FAIL_EXPQI_FLOW')   # every information set routes integrally: s >= Wm for every frame (C1)
