#!/usr/bin/env python3
"""expRJ (E2): cut_ceiling.py extended to the repo's charge.

Ceiling (zll-recursion dictionary, expQC §2; expL7 §3D):
    rho = s/(Wm) >= (routing rate)/(coding rate) >= lambda*/sigma
    a_cap = ln(sigma_hat/lambda)/ln b
  - lambda* = max CONCURRENT multicommodity-flow rate (unit capacities = one symbol per physical edge).
    Exact LP (HiGHS) for <= 512 roles. Otherwise lambda_LB = 1/maxload of an explicit ECMP (shortest-path split)
    flow; that is a feasible flow, so a_cap(ECMP) >= a_cap(LP): still a valid ceiling.
    lambda_vol = C/sum dist (volume bound) is an UPPER bound on lambda*, so a_vol is a lower end of the ceiling.
  - sigma_hat = best cut found (cut capacity / demand crossing): an upper bound on the undirected coding rate.
  - b = own-pair distance; demands are kept only if dist >= b0 and b = b0 (uniform metric x = 1/b0 is then a valid dual).
  - expL7's uniform-hop form a_hop = ln(k b sigma/E)/ln b uses lambda ~ E/(kb) >= lambda* (it under-states the ceiling).
Rank charge (t-ary MDS square gates): every port is one physical symbol; a far-pair role at rate r occupies r terminal
ports at its gate, so arity t = d_int + r (fractional r allowed as a generous symbol-extension charge); a Tanner session
is one graph edge (one port at each end), arity t = deg. Non-square [t,k!=t/2] MDS gates are not bijections and are not
admissible in the contract; for square gates Sum_g (t - dim) = E, so the expQK S/E charge reduces to edges.
"""
import math, sys, time, json
import numpy as np
import numba as nb
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

# ---------------------------------------------------------------- graphs (edge lists u<v, cap)
def _finish(n, edges):
    d = {}
    for u, v in edges:
        if u == v: continue
        if u > v: u, v = v, u
        d[(u, v)] = d.get((u, v), 0) + 1
    E = np.array(sorted(d), dtype=np.int64).reshape(-1, 2)
    cap = np.array([d[tuple(e)] for e in E], dtype=np.float64)
    return n, E, cap

def rand_regular(n, d, seed):
    import networkx as nx
    while True:
        G = nx.random_regular_graph(d, n, seed=seed)
        if nx.is_connected(G): break
        seed += 1000
    return _finish(n, list(G.edges()))

def rand_biregular(n1, t, D, seed):
    """bipartite physical graph: n1 gates of arity t, n2 = n1 t / D gates of arity D (simple, by swaps)."""
    rng = np.random.default_rng(seed)
    n2 = n1 * t // D; assert n1 * t == n2 * D
    L = np.repeat(np.arange(n1), t); R = np.repeat(np.arange(n2), D); rng.shuffle(R)
    # repair parallel edges by random swaps (multi-edges that survive are merged into capacity)
    for _ in range(50):
        seen = {}; dup = []
        for j, (u, v) in enumerate(zip(L, R)):
            if (u, v) in seen: dup.append(j)
            else: seen[(u, v)] = j
        if not dup: break
        for j in dup:
            o = rng.integers(len(R)); R[j], R[o] = R[o], R[j]
    return _finish(n1 + n2, [(int(u), n1 + int(v)) for u, v in zip(L, R)])

def lps(p, q):
    """LPS X^{p,q}: Cayley graph of PSL/PGL(2,q), (p+1)-regular. p,q primes = 1 mod 4."""
    i = next(x for x in range(q) if (x * x + 1) % q == 0)
    sols = []
    r = int(math.isqrt(p))
    for a0 in range(1, r + 1, 2):
        for a1 in range(-r, r + 1):
            for a2 in range(-r, r + 1):
                for a3 in range(-r, r + 1):
                    if a1 % 2 or a2 % 2 or a3 % 2: continue
                    if a0 * a0 + a1 * a1 + a2 * a2 + a3 * a3 == p:
                        sols.append(((a0 + i * a1) % q, (a2 + i * a3) % q, (-a2 + i * a3) % q, (a0 - i * a1) % q))
    assert len(sols) == p + 1, len(sols)
    def norm(m):
        a, b, c, d = m
        f = a if a else b
        inv = pow(f, q - 2, q)
        return (a * inv % q, b * inv % q, c * inv % q, d * inv % q)
    def mul(m, g):
        a, b, c, d = m; e, f, h, k = g
        return norm(((a * e + b * h) % q, (a * f + b * k) % q, (c * e + d * h) % q, (c * f + d * k) % q))
    idx = {norm((1, 0, 0, 1)): 0}; order = [norm((1, 0, 0, 1))]; edges = []
    h = 0
    while h < len(order):
        m = order[h]
        for g in sols:
            w = mul(m, g)
            if w not in idx:
                idx[w] = len(order); order.append(w)
            if h < idx[w]: edges.append((h, idx[w]))
        h += 1
    return _finish(len(order), edges)

def margulis(n):
    """Gabber-Galil / Margulis on Z_n^2: 8-regular multigraph (loops dropped)."""
    V = lambda x, y: (x % n) * n + (y % n)
    edges = []
    for x in range(n):
        for y in range(n):
            u = V(x, y)
            for w in (V(x + 2 * y, y), V(x + 2 * y + 1, y), V(x, y + 2 * x), V(x, y + 2 * x + 1)):
                edges.append((u, w))
    return _finish(n * n, edges)

def csr(n, E, cap):
    m = len(E)
    src = np.concatenate([E[:, 0], E[:, 1]]); dst = np.concatenate([E[:, 1], E[:, 0]])
    eid = np.concatenate([np.arange(m), np.arange(m)])
    o = np.lexsort((dst, src))
    src, dst, eid = src[o], dst[o], eid[o]
    indptr = np.zeros(n + 1, np.int64); np.add.at(indptr, src + 1, 1); indptr = np.cumsum(indptr)
    return indptr, dst.astype(np.int64), eid.astype(np.int64)

def dcsr(n, P):
    """demand graph CSR (a vertex may carry several demands)."""
    Ed = np.asarray(P[:, :2], dtype=np.int64)
    ip, ind, _ = csr(n, Ed, np.ones(len(Ed)))
    return ip, ind

# ---------------------------------------------------------------- numba kernels
@nb.njit(cache=True)
def bfs(indptr, ind, s, dist, order):
    n = len(dist)
    for i in range(n): dist[i] = -1
    dist[s] = 0; order[0] = s; h = 0; t = 1
    while h < t:
        u = order[h]; h += 1
        for j in range(indptr[u], indptr[u + 1]):
            w = ind[j]
            if dist[w] < 0:
                dist[w] = dist[u] + 1; order[t] = w; t += 1
    return t

@nb.njit(cache=True)
def far_pairs(indptr, ind, n, b0, seed):
    """greedy: random order, pair u with a random free v at dist >= b0 (prefer exactly the max available >= b0)."""
    np.random.seed(seed)
    perm = np.random.permutation(n)
    free = np.ones(n, np.bool_)
    dist = np.empty(n, np.int64); order = np.empty(n, np.int64)
    P = np.empty((n // 2, 3), np.int64); k = 0
    cand = np.empty(n, np.int64)
    for u in perm:
        if not free[u]: continue
        t = bfs(indptr, ind, u, dist, order)
        c = 0
        for j in range(t):
            w = order[j]
            if free[w] and w != u and dist[w] >= b0:
                cand[c] = w; c += 1
        if c == 0: continue
        v = cand[np.random.randint(c)]
        free[u] = False; free[v] = False
        P[k, 0] = u; P[k, 1] = v; P[k, 2] = dist[v]; k += 1
    return P[:k]

@nb.njit(cache=True)
def ecmp_load(indptr, ind, eid, n, m, P):
    load = np.zeros(m)
    dist = np.empty(n, np.int64); order = np.empty(n, np.int64)
    sig = np.zeros(n); fl = np.zeros(n)
    for q in range(P.shape[0]):
        s = P[q, 0]; tt = P[q, 1]
        t = bfs(indptr, ind, s, dist, order)
        for j in range(t): sig[order[j]] = 0.0
        sig[s] = 1.0
        for j in range(t):
            u = order[j]
            for a in range(indptr[u], indptr[u + 1]):
                w = ind[a]
                if dist[w] == dist[u] + 1: sig[w] += sig[u]
        for j in range(t): fl[order[j]] = 0.0
        fl[tt] = 1.0
        for j in range(t - 1, -1, -1):
            v = order[j]
            if fl[v] == 0.0 or v == s: continue
            for a in range(indptr[v], indptr[v + 1]):
                u = ind[a]
                if dist[u] == dist[v] - 1:
                    f = fl[v] * sig[u] / sig[v]
                    load[eid[a]] += f; fl[u] += f
    return load

@nb.njit(cache=True)
def local_search(indptr, ind, eid, cap, dptr, dind, side, passes):
    n = len(side)
    cut = 0.0; cross = 0
    for u in range(n):
        for a in range(indptr[u], indptr[u + 1]):
            if side[ind[a]] != side[u]: cut += cap[eid[a]]
        for a in range(dptr[u], dptr[u + 1]):
            if side[dind[a]] != side[u]: cross += 1
    cut /= 2; cross //= 2
    best = cut / cross if cross > 0 else 1e18
    for _ in range(passes):
        improved = False
        for u in range(n):
            dc = 0.0
            for a in range(indptr[u], indptr[u + 1]):
                dc += cap[eid[a]] if side[ind[a]] == side[u] else -cap[eid[a]]
            dx = 0
            for a in range(dptr[u], dptr[u + 1]):
                dx += 1 if side[dind[a]] == side[u] else -1
            nc = cut + dc; nx = cross + dx
            if nx > 0 and nc / nx < best - 1e-12:
                side[u] = 1 - side[u]; cut = nc; cross = nx; best = nc / nx; improved = True
        if not improved: break
    return best

@nb.njit(cache=True)
def sweep(indptr, ind, eid, cap, dptr, dind, order):
    n = len(order)
    side = np.zeros(n, np.int64)
    cut = 0.0; cross = 0; best = 1e18; bi = -1
    for i in range(n - 1):
        v = order[i]; side[v] = 1
        for a in range(indptr[v], indptr[v + 1]):
            cut += -cap[eid[a]] if side[ind[a]] == 1 else cap[eid[a]]
        for a in range(dptr[v], dptr[v + 1]):
            cross += -1 if side[dind[a]] == 1 else 1
        if cross > 0 and cut / cross < best:
            best = cut / cross; bi = i
    return best, bi

# ---------------------------------------------------------------- sparsest-cut upper bound
def sigma_hat(n, E, cap, P, seed=0, nvec=4, nballs=24, passes=30):
    indptr, ind, eid = csr(n, E, cap)
    dptr, dind = dcsr(n, P)
    deg = np.zeros(n); np.add.at(deg, E[:, 0], cap); np.add.at(deg, E[:, 1], cap)
    nd = np.diff(dptr); best = float((deg[nd > 0] / nd[nd > 0]).min())  # singletons
    A = sp.coo_matrix((np.concatenate([cap, cap]), (np.concatenate([E[:, 0], E[:, 1]]), np.concatenate([E[:, 1], E[:, 0]]))), shape=(n, n)).tocsr()
    c = 2 * deg.max()
    M = sp.diags(c - deg) + A
    rng = np.random.default_rng(seed)
    try:
        vals, vecs = eigsh(M, k=nvec + 1, which='LA', tol=1e-6, maxiter=20000, v0=rng.standard_normal(n))
        vecs = vecs[:, np.argsort(-vals)][:, 1:]
    except Exception:
        vecs = rng.standard_normal((n, nvec))
    dirs = [vecs[:, j] for j in range(vecs.shape[1])]
    for _ in range(6):
        dirs.append(vecs @ rng.standard_normal(vecs.shape[1]))
    dist = np.empty(n, np.int64); order = np.empty(n, np.int64)
    for _ in range(nballs):
        s = rng.integers(n); bfs(indptr, ind, s, dist, order)
        dirs.append(dist + 0.01 * rng.random(n))
    for f in dirs:
        o = np.argsort(f, kind='stable').astype(np.int64)
        val, bi = sweep(indptr, ind, eid, cap, dptr, dind, o)
        if bi >= 0:
            side = np.zeros(n, np.int64); side[o[:bi + 1]] = 1
            val = min(val, local_search(indptr, ind, eid, cap, dptr, dind, side, passes))
        best = min(best, val)
    return best

# ---------------------------------------------------------------- exact max concurrent flow (HiGHS)
def mcf_lp(n, E, cap, P, time_limit=1500):
    from scipy.optimize import linprog
    m = len(E); k = len(P); A2 = 2 * m
    tail = np.concatenate([E[:, 0], E[:, 1]]); head = np.concatenate([E[:, 1], E[:, 0]])
    nv = k * A2 + 1
    # conservation rows (i, v): sum_out - sum_in - lam*(v==s) + lam*(v==t) = 0 ; drop row v = t_i (redundant)
    rows = []; cols = []; vals = []
    for i in range(k):
        base = i * A2; rb = i * n
        a = np.arange(A2)
        rows += [rb + tail, rb + head]; cols += [base + a, base + a]; vals += [np.ones(A2), -np.ones(A2)]
        rows.append(np.array([rb + P[i, 0]])); cols.append(np.array([nv - 1])); vals.append(np.array([-1.0]))
    rows = np.concatenate(rows); cols = np.concatenate(cols); vals = np.concatenate(vals)
    keep = np.ones(k * n, bool); keep[np.arange(k) * n + P[:, 1]] = False
    remap = -np.ones(k * n, np.int64); remap[keep] = np.arange(keep.sum())
    msk = keep[rows]
    Aeq = sp.csr_matrix((vals[msk], (remap[rows[msk]], cols[msk])), shape=(keep.sum(), nv))
    beq = np.zeros(keep.sum())
    # capacity rows
    ce = np.concatenate([np.arange(m), np.arange(m)])
    r2 = np.tile(ce, k); c2 = (np.arange(k)[:, None] * A2 + np.arange(A2)[None, :]).ravel()
    Aub = sp.csr_matrix((np.ones(len(r2)), (r2, c2)), shape=(m, nv))
    cvec = np.zeros(nv); cvec[-1] = -1.0
    t0 = time.time()
    res = linprog(cvec, A_ub=Aub, b_ub=cap, A_eq=Aeq, b_eq=beq, bounds=(0, None), method='highs-ipm',
                  options=dict(time_limit=time_limit, presolve=True))
    lam = -res.fun if res.status == 0 else float('nan')
    dual = -np.asarray(res.ineqlin.marginals) if res.status == 0 else None
    return lam, time.time() - t0, res.status, dual

def dual_cuts(n, E, cap, P, ell):
    """Leighton-Rao style rounding: sweep balls of the LP dual metric around every terminal."""
    from scipy.sparse.csgraph import dijkstra
    indptr, ind, eid = csr(n, E, cap)
    dptr, dind = dcsr(n, P)
    w = np.maximum(ell, 0) + 1e-9
    A = sp.coo_matrix((np.concatenate([w, w]), (np.concatenate([E[:, 0], E[:, 1]]), np.concatenate([E[:, 1], E[:, 0]]))), shape=(n, n)).tocsr()
    src = np.unique(P[:, :2].ravel())
    D = dijkstra(A, indices=src)
    best = 1e18
    for row in D:
        o = np.argsort(row, kind='stable').astype(np.int64)
        val, bi = sweep(indptr, ind, eid, cap, dptr, dind, o)
        if bi >= 0:
            side = np.zeros(n, np.int64); side[o[:bi + 1]] = 1
            val = min(val, local_search(indptr, ind, eid, cap, dptr, dind, side, 30))
        best = min(best, val)
    return best

@nb.njit(cache=True)
def exact_sigma(indptr, ind, eid, cap, dptr, dind):
    """brute force over all 2^(n-1) cuts (Gray code), n <= ~26."""
    n = len(dptr) - 1
    side = np.zeros(n, np.int64); cut = 0.0; cross = 0; best = 1e18
    N = 1 << (n - 1)
    for g in range(1, N):
        v = 0
        x = g
        while (x & 1) == 0:
            x >>= 1; v += 1
        for a in range(indptr[v], indptr[v + 1]):
            cut += -cap[eid[a]] if side[ind[a]] != side[v] else cap[eid[a]]
        for a in range(dptr[v], dptr[v + 1]):
            cross += -1 if side[dind[a]] != side[v] else 1
        side[v] = 1 - side[v]
        if cross > 0:
            r = cut / cross
            if r < best: best = r
    return best

# ---------------------------------------------------------------- one instance
def evaluate(name, n, E, cap, P, b0, rate_cap=None, arity_base=None, lp=False, seed=0):
    m = len(E); k = len(P)
    indptr, ind, eid = csr(n, E, cap)
    C = cap.sum()
    sig = sigma_hat(n, E, cap, P, seed=seed)
    load = ecmp_load(indptr, ind, eid, n, m, P)
    lam_ecmp = 1.0 / (load / cap).max()
    lam_vol = C / P[:, 2].sum()
    lam_lp = float('nan'); tlp = 0.0
    sig_heur = sig
    if lp:
        lam_lp, tlp, st, dual = mcf_lp(n, E, cap, P)
        if dual is not None:
            sig = min(sig, dual_cuts(n, E, cap, P, dual))
    if n <= 26:
        dptr, dind = dcsr(n, P)
        sig = exact_sigma(indptr, ind, eid, cap, dptr, dind)
    r = sig if rate_cap is None else min(sig, rate_cap)
    lam = lam_lp if lp and lam_lp == lam_lp else lam_ecmp
    b = b0; bmean = P[:, 2].mean()
    def A(l, bb): return math.log(r / l) / math.log(bb)
    out = dict(name=name, n=n, E=float(C), k=k, b0=b0, bmean=round(float(bmean), 3), sigma=round(sig, 4), sigma_heur=round(sig_heur, 4), exact_sigma=(n <= 26), rate=round(r, 4),
               lam_vol=round(lam_vol, 5), lam_ecmp=round(lam_ecmp, 5), lam_lp=round(lam_lp, 5) if lp else None,
               gap=round(r / lam, 3), a_cap=round(A(lam, b), 4), a_vol=round(A(lam_vol, b), 4),
               a_hop_expL7=round(math.log(k * bmean * r / C) / math.log(bmean), 4),
               a_cap_bmean=round(A(lam, bmean), 4), tlp=round(tlp, 1))
    if arity_base is not None:
        out['arity'] = round(arity_base(r), 2)
    return out

@nb.njit(cache=True)
def far_pairs_fast(indptr, ind, n, b0, seed, tries):
    """pair u with a random free vertex outside its (b0-1)-ball (exact dist >= b0, recorded by a full BFS)."""
    np.random.seed(seed)
    perm = np.random.permutation(n)
    free = np.ones(n, np.bool_); freelist = perm.copy(); nfree = n
    pos = np.empty(n, np.int64)
    for i in range(n): pos[freelist[i]] = i
    dist = np.empty(n, np.int64); order = np.empty(n, np.int64)
    P = np.empty((n // 2, 3), np.int64); k = 0
    mark = -np.ones(n, np.int64)
    for u in perm:
        if not free[u]: continue
        # ball of radius b0-1
        mark[u] = u; order[0] = u; dist[u] = 0; h = 0; t = 1
        while h < t:
            x = order[h]; h += 1
            if dist[x] >= b0 - 1: continue
            for j in range(indptr[x], indptr[x + 1]):
                w = ind[j]
                if mark[w] != u:
                    mark[w] = u; dist[w] = dist[x] + 1; order[t] = w; t += 1
        v = -1
        for _ in range(tries):
            if nfree <= 1: break
            c = freelist[np.random.randint(nfree)]
            if c != u and mark[c] != u:
                v = c; break
        if v < 0: continue
        for z in (u, v):
            free[z] = False; i = pos[z]; last = freelist[nfree - 1]
            freelist[i] = last; pos[last] = i; nfree -= 1
        P[k, 0] = u; P[k, 1] = v; P[k, 2] = b0; k += 1
    return P[:k]

def a_H(r, k, b0, E):
    """uniform-hop dictionary ceiling, consistent form: a = 1 - ln(E/(r k))/ln b_min."""
    return math.log(r * k * b0 / E) / math.log(b0)
