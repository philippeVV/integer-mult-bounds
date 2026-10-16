"""Bootstrap upper bounds by simulated annealing over vertex orders (float-key representation).
cost(sigma) = sum_v (k_v - j_v)^+, j_v = #neighbours earlier.  Boot = min over orders (any order = valid upper bound).
Also: spectral lower bound Boot >= lambda2 * n / 8 (bisection at the middle prefix; cost >= max prefix cut / 2... see README).
usage: bootsa.py GRAPHFILE NAME restarts sweeps [seed]   -> JSON line"""
import sys, json, numpy as np, numba as nb, time

def load(gf, name):
    G = json.load(open(gf))[name]
    V = G['Pts'] + G['Lns']; E = G['E']
    ends = [[] for _ in range(E)]
    for v, es in enumerate(V):
        for e in es: ends[e].append(v)
    n = len(V)
    deg = np.array([len(es) for es in V], dtype=np.int64)
    ptr = np.zeros(n + 1, dtype=np.int64); ptr[1:] = np.cumsum(deg)
    adj = np.zeros(ptr[-1], dtype=np.int64); fill = ptr[:-1].copy()
    for a, b in ends:
        adj[fill[a]] = b; fill[a] += 1; adj[fill[b]] = a; fill[b] += 1
    return G, n, E, ptr, adj, deg

@nb.njit(cache=True)
def cost_of(key, ptr, adj, k):
    n = len(key); c = 0; j = np.zeros(n, np.int64)
    for v in range(n):
        for t in range(ptr[v], ptr[v + 1]):
            if key[adj[t]] < key[v]: j[v] += 1
        c += max(k[v] - j[v], 0)
    return c, j

@nb.njit(cache=True)
def anneal(key, ptr, adj, k, sweeps, T0, T1, seed):
    np.random.seed(seed)
    n = len(key)
    c, j = cost_of(key, ptr, adj, k)
    best = c; bestkey = key.copy()
    nsteps = sweeps * n
    for it in range(nsteps):
        T = T0 * (T1 / T0) ** (it / nsteps)
        v = np.random.randint(n)
        r = np.random.random()
        if r < 0.5:
            # move next to a random neighbour-of-neighbour or neighbour
            u = adj[ptr[v] + np.random.randint(ptr[v + 1] - ptr[v])]
            if np.random.random() < 0.5:
                u = adj[ptr[u] + np.random.randint(ptr[u + 1] - ptr[u])]
            newk = key[u] + (np.random.random() - 0.5) * 1e-9
        elif r < 0.8:
            w = np.random.randint(n)
            newk = key[w] + (np.random.random() - 0.5) * 1e-9
        else:
            newk = np.random.random()
        old = key[v]
        lo = min(old, newk); hi = max(old, newk)
        # neighbours u in (lo, hi)
        dj_v = 0; delta = 0
        for t in range(ptr[v], ptr[v + 1]):
            u = adj[t]
            if key[u] > lo and key[u] < hi:
                if newk > old:  # v moves later past u: j_v+1, j_u-1
                    dj_v += 1
                    delta += max(k[u] - (j[u] - 1), 0) - max(k[u] - j[u], 0)
                else:
                    dj_v -= 1
                    delta += max(k[u] - (j[u] + 1), 0) - max(k[u] - j[u], 0)
        delta += max(k[v] - (j[v] + dj_v), 0) - max(k[v] - j[v], 0)
        if delta <= 0 or np.random.random() < np.exp(-delta / T):
            for t in range(ptr[v], ptr[v + 1]):
                u = adj[t]
                if key[u] > lo and key[u] < hi:
                    if newk > old: j[u] -= 1
                    else: j[u] += 1
            j[v] += dj_v; key[v] = newk; c += delta
            if c < best:
                best = c; bestkey[:] = key
    return best, bestkey

def renorm(key):
    o = np.argsort(key, kind='stable'); out = np.empty_like(key); out[o] = (np.arange(len(key)) + 0.5) / len(key); return out

def greedy_key(n, ptr, adj, k, rng):
    j = np.zeros(n, np.int64); placed = np.zeros(n, bool); key = np.zeros(n)
    for pos in range(n):
        rem = np.where(placed, 10 ** 9, np.maximum(k - j, 0))
        m = rem.min(); cand = np.nonzero(rem == m)[0]; v = rng.choice(cand)
        placed[v] = True; key[v] = (pos + 0.5) / n
        for t in range(ptr[v], ptr[v + 1]): j[adj[t]] += 1
    return key

def spectral_lb(n, ptr, adj, deg):
    import scipy.sparse as sp, scipy.sparse.linalg as sla
    rows = np.repeat(np.arange(n), deg); A = sp.csr_matrix((np.ones(len(adj)), (rows, adj)), shape=(n, n))
    L = sp.diags(deg.astype(float)) - A
    if n <= 3000:
        w = np.linalg.eigvalsh(L.toarray()); lam2 = w[1]
    else:
        w = sla.eigsh(L, k=3, sigma=-1e-3, which='LM', return_eigenvectors=False); lam2 = sorted(w)[1]
    # any order: cost >= (1/2) max prefix cut >= (1/2) cut(prefix of size floor(n/2)) >= lam2 * m(n-m)/(2n)
    m = n // 2
    return float(lam2), lam2 * m * (n - m) / (2 * n)

if __name__ == '__main__':
    gf, name, R, sweeps = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    G, n, E, ptr, adj, deg = load(gf, name)
    k = deg // 2
    rng = np.random.default_rng(seed)
    t0 = time.time()
    res = {}
    # greedy
    gbest = 10 ** 9; gkey = None
    for r in range(max(3, R)):
        key = greedy_key(n, ptr, adj, k, rng); c, _ = cost_of(key, ptr, adj, k)
        if c < gbest: gbest, gkey = c, key
    res['greedy'] = int(gbest)
    best = gbest; bkey = gkey
    # cyclic orders for Singer graphs: point i at i, line j at j + shift
    nP = len(G['Pts'])
    if nP == n - nP:
        cyc = 10 ** 9
        for sh in np.linspace(0, nP, 4 * nP, endpoint=False):
            key = np.concatenate([np.arange(nP), np.arange(nP) + sh + 1e-3]).astype(float)
            key = renorm(key); c, _ = cost_of(key, ptr, adj, k)
            if c < cyc: cyc, ckey = c, key
        res['cyclic'] = int(cyc)
        if cyc < best: best, bkey = cyc, ckey
    for r in range(R):
        start = bkey.copy() if r % 2 == 0 else renorm(rng.random(n))
        T0 = 0.6 if r % 2 == 0 else 2.0
        c, key = anneal(start, ptr, adj, k, sweeps, T0, 0.05, seed * 1000 + r)
        key = renorm(key); c2, _ = cost_of(key, ptr, adj, k); assert c2 == c, (c, c2)
        if c < best: best, bkey = c, key
    lam2, lb = spectral_lb(n, ptr, adj, deg)
    res.update(graph=name, n=n, E=E, d=[int(deg.min()), int(deg.max())], girth=G.get('girth'), boot_ub=int(best), ub_over_E=round(best / E, 5),
               lam2=round(lam2, 4), spec_lb=round(lb, 2), lb_over_E=round(lb / E, 5), sec=round(time.time() - t0, 1))
    np.save(f'bkey_{name}.npy', bkey)
    print(json.dumps(res), flush=True)
