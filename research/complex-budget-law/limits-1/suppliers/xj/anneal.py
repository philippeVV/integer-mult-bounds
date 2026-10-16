"""Simulated annealing upper bound for min s (phase-frame beta model) on BH scenario A.
usage: anneal.py m {cost|rank} restarts sweeps seed"""
import sys, os, random, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
import beta as B

m, var, R, SW, seed = int(sys.argv[1]), sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
nb, rk, cost, I = B.tables(m)
tab = cost if var == 'cost' else rk
nbr = [[] for _ in range(B.NV)]           # (other, offset)
for u, v in B.PHYS: nbr[u].append((v, 0)); nbr[v].append((u, 0))
for a, b in B.SESS: nbr[a].append((b, I)); nbr[b].append((a, I))
rng = random.Random(seed); NM = 1 << nb
best_all = None
for r in range(R):
    beta = [rng.randrange(NM) if r % 2 else 0 for _ in range(B.NV)]
    s = B.s_of(beta, tab, I); best = s; bb = beta[:]
    T0, T1 = 2.0, 0.05
    steps = SW * B.NV
    for it in range(steps):
        T = T0 * (T1 / T0) ** (it / steps)
        v = rng.randrange(B.NV); new = rng.randrange(NM)
        if new == beta[v]: continue
        old = beta[v]; d = 0
        for w, off in nbr[v]:
            d += tab[new ^ beta[w] ^ off] - tab[old ^ beta[w] ^ off]
        if d <= 0 or rng.random() < math.exp(-d / T):
            beta[v] = new; s += d
            if s < best: best, bb = s, beta[:]
    assert B.s_of(bb, tab, I) == best
    print(f"m={m} {var} restart {r}: best s = {best}  (Wm = {157*m}, flow {147*m})", flush=True)
    if best_all is None or best < best_all[0]: best_all = (best, bb)
print(f"FINAL m={m} {var}: min found s = {best_all[0]}  Wm = {157*m}  delta = {157*m - best_all[0]}")
print("beta =", best_all[1])
