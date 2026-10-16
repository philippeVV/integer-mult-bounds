#!/usr/bin/env python3
"""expRJ Tanner-hybrid table: girth of LPS/Margulis/PG/random seeds and a_T = 1 - ln((1-s)/s)/ln(g-1)
at s = S/E in {1/4 (Hermitian asymptote, expQK), 1/2 (cut bound for balanced sessions)}."""
import sys, math, json
sys.path.insert(0, '.')
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from cc2 import *
from sweep import pg2

@nb.njit(cache=True)
def girth_from(indptr, ind, s, n):
    dist = -np.ones(n, np.int64); par = -np.ones(n, np.int64); q = np.empty(n, np.int64)
    dist[s] = 0; q[0] = s; h = 0; t = 1; best = 1 << 30
    while h < t:
        u = q[h]; h += 1
        if 2 * dist[u] + 1 >= best: break
        for j in range(indptr[u], indptr[u + 1]):
            w = ind[j]
            if dist[w] < 0:
                dist[w] = dist[u] + 1; par[w] = u; q[t] = w; t += 1
            elif w != par[u]:
                best = min(best, dist[u] + dist[w] + 1)
    return best

def girth(G, samples):
    n, E, cap = G
    if (cap > 1).any(): return 2
    ip, ind, eid = csr(n, E, cap)
    rng = np.random.default_rng(0)
    src = range(n) if samples is None else rng.integers(0, n, samples)
    return min(girth_from(ip, ind, int(s), n) for s in src)

def aT(s, g): return 1 - math.log((1 - s) / s) / math.log(g - 1)
rows = []
for name, G, deg, samp in [("PG(2,5) inc", pg2(5), 6, None), ("PG(2,13) inc", pg2(13), 14, 50),
                           ("LPS(13,5)", lps(13, 5), 14, 1), ("LPS(5,13)", lps(5, 13), 6, 1), ("LPS(5,17)", lps(5, 17), 6, 1),
                           ("LPS(5,29)", lps(5, 29), 6, 1), ("LPS(5,37)", lps(5, 37), 6, 1), ("LPS(13,17)", lps(13, 17), 14, 1),
                           ("Margulis64", margulis(64), 8, 200), ("rr4 n=4096", rand_regular(4096, 4, 1), 4, None),
                           ("rr6 n=4096", rand_regular(4096, 6, 1), 6, None)]:
    n, E, cap = G
    g = girth(G, samp)
    Etot = float(cap.sum())
    r = dict(seed=name, n=n, deg=deg, E=Etot, girth=g, b=g - 1, S_at_quarter=round(Etot / 4),
             a_T_quarter=round(aT(0.25, g), 4) if g > 2 else None, a_T_half_cutbound=1.0 if g > 2 else None)
    print(json.dumps(r), flush=True); rows.append(r)
print("analytic a_T(1/4, g):", {g: round(aT(0.25, g), 4) for g in (6, 8, 10, 12, 14, 16)})
print("girth needed for a>=1/2 at S/E=s:", {s: math.ceil(1 + ((1 - s) / s) ** 2) for s in (0.177, 0.2, 0.25, 0.3)})
json.dump(rows, open('girth.json', 'w'), indent=1)
