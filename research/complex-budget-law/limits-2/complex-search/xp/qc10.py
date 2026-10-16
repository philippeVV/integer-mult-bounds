"""Quasi-cyclic (voltage) lifts of K_{d,d} with girth >= G: exponent matrix B (d x d) over Z_L, check (i,x) ~ var (j, x+B_ij).
Randomised greedy fill with girth check by BFS.  usage: qc10.py d Lmin Lmax tries [G=10]  -> qc10.json entries QC{d}_{L}"""
import sys, json, numpy as np, numba as nb, os

@nb.njit(cache=True)
def girth_from(src, ptr, adj, n, cap):
    dist = -np.ones(n, np.int64); par = -np.ones(n, np.int64)
    q = np.empty(n, np.int64); h = 0; t = 0; q[t] = src; t += 1; dist[src] = 0; best = cap
    while h < t:
        u = q[h]; h += 1
        if 2 * dist[u] + 1 >= best: break
        for a in range(ptr[u], ptr[u + 1]):
            w = adj[a]
            if w == par[u]: continue
            if dist[w] < 0:
                dist[w] = dist[u] + 1; par[w] = u; q[t] = w; t += 1
            else:
                c = dist[u] + dist[w] + 1
                if c < best: best = c
    return best

def build(B, mask, d, L):
    # vertices: checks i*L+x (0..dL-1), vars dL + j*L + y
    n = 2 * d * L; nbrs = [[] for _ in range(n)]
    for i in range(d):
        for j in range(d):
            if not mask[i, j]: continue
            for x in range(L):
                c = i * L + x; v = d * L + j * L + (x + B[i, j]) % L
                nbrs[c].append(v); nbrs[v].append(c)
    deg = np.array([len(a) for a in nbrs]); ptr = np.zeros(n + 1, np.int64); ptr[1:] = np.cumsum(deg)
    adj = np.array([w for a in nbrs for w in a], dtype=np.int64)
    return n, ptr, adj

def girth(B, mask, d, L, cap=99):
    n, ptr, adj = build(B, mask, d, L); g = cap
    for i in range(d):  # by circulant symmetry, one copy per base vertex suffices
        g = min(g, girth_from(i * L, ptr, adj, n, cap)); g = min(g, girth_from(d * L + i * L, ptr, adj, n, cap))
    return g

def attempt(d, L, G, rng, per_entry=200):
    B = np.zeros((d, d), np.int64); mask = np.zeros((d, d), bool)
    mask[0, :] = True; mask[:, 0] = True  # first row/col zero WLOG
    order = [(i, j) for i in range(1, d) for j in range(1, d)]
    for (i, j) in order:
        mask[i, j] = True; ok = False
        for _ in range(per_entry):
            B[i, j] = rng.integers(1, L)
            if girth(B, mask, d, L, G + 2) >= G: ok = True; break
        if not ok: return None
    return B

if __name__ == '__main__':
    d, L0, L1, tries = map(int, sys.argv[1:5]); G = int(sys.argv[5]) if len(sys.argv) > 5 else 10
    rng = np.random.default_rng(1); out = json.load(open('qc10.json')) if os.path.exists('qc10.json') else {}
    step = int(os.environ.get("STEP", "1"))
    for L in range(L0, L1 + 1, step):
        found = None
        for t in range(tries):
            B = attempt(d, L, G, rng)
            if B is not None: found = B; break
        if found is None:
            print(json.dumps(dict(d=d, L=L, found=False)), flush=True); continue
        g = girth(found, np.ones((d, d), bool), d, L)
        # write graph: Pts = checks, Lns = vars, edge ids
        Pts = [[] for _ in range(d * L)]; Lns = [[] for _ in range(d * L)]; E = 0
        for i in range(d):
            for j in range(d):
                for x in range(L):
                    Pts[i * L + x].append(E); Lns[j * L + (x + found[i, j]) % L].append(E); E += 1
        out[f'QC{d}_{L}'] = dict(Pts=Pts, Lns=Lns, E=E, girth=int(g), B=found.tolist(), L=L)
        json.dump(out, open('qc10.json', 'w'))
        print(json.dumps(dict(d=d, L=L, found=True, girth=int(g), E=E, n=2 * d * L, B=found.tolist())), flush=True)
        break
