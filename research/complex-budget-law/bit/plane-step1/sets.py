"""Enumerate maximal pairwise-odd common-r sets as in probe_plane.py; keep (|S|, rankM) classes."""
import core, json
from collections import Counter
n = 111; live = core.LINES[:111]
lines = list(range(120))
out = []; stats = Counter()
for r in lines:
    perp = [a for a in range(n) if core.G[a][r] == 0 and a != r]
    T = {a: sum(1 << j for j in range(n) if j != a and core.G[a][j] == 0 and core.G[j][r] != 0) for a in perp}
    adj = {a: {b for b in perp if b != a and core.G[a][b] != 0} for a in perp}
    cl = []
    def bk(Rs, P, X):
        if not P and not X: cl.append(sorted(Rs)); return
        for v in list(P):
            bk(Rs | {v}, P & adj[v], X & adj[v]); P = P - {v}; X = X | {v}
    bk(set(), set(perp), set())
    for S in cl:
        rk = core.rank2([T[a] for a in S])
        stats[(len(S), 2 * rk)] += 1
        out.append((r, S, rk))
print(sorted(stats.items()))
json.dump(out, open('sets_all.json', 'w'))
