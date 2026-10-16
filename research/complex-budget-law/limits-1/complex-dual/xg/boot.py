"""Gauge-degeneration (= bootstrap) bound, any field: dim W <= min over vertex orders sum_v (deg(v)/2 - #earlier nbrs)^+.
Proof in README sec 2 (limit of the gauge orbit t^c, upper semicontinuity of dim of intersections)."""
import random, json, sys
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from xg_core import load
def boot(V, E, R=300, seed=1):
    ends = [[] for _ in range(E)]
    for v, es in enumerate(V):
        for e in es: ends[e].append(v)
    n = len(V); adj = [[] for _ in range(n)]
    for a, b in ends: adj[a].append(b); adj[b].append(a)
    k = [len(es) // 2 for es in V]; rng = random.Random(seed); best = None
    for t in range(R):
        j = [0] * n; placed = [False] * n; c = 0; v = rng.randrange(n); order = []
        while True:
            placed[v] = True; order.append(v); c += max(k[v] - j[v], 0)
            for w in adj[v]: j[w] += 1
            if len(order) == n: break
            bestd = min(max(k[w] - j[w], 0) for w in range(n) if not placed[w])
            v = rng.choice([w for w in range(n) if not placed[w] and max(k[w] - j[w], 0) == bestd])
        if best is None or c < best[0]: best = (c, order)
    return best
if __name__ == '__main__':
    for name in sys.argv[1:]:
        V, nP, E, g = load(name); c, o = boot(V, E)
        print(json.dumps(dict(graph=name, E=E, girth=g, deg=len(V[0]), boot=c, ratio=round(c / E, 4), threshold=round(1 / g, 4))))
