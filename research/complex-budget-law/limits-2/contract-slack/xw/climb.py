"""Hill-climb maximising greedy-Boot over orientations, started from the most unbalanced class profile
(t=1: repetition points, parity lines). Walk accepts -1 moves w.p. 0.3. Tests whether the balanced (half-rate) profile is the attractor."""
import random, sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # python -I omits the script dir
from xw import pg, boot, orient
for q in (3, 4, 5):
    N, edges, adj = pg(q); E = len(edges); d = q + 1; n = 2 * N
    k = [1] * N + [q] * N; arcs = orient(edges, n, k); cur = boot(adj, k, 20, 11); hr = random.Random(7)
    start = cur; best = (cur, list(k))
    for it in range(4000):
        i = hr.randrange(E); a, b = arcs[i]
        k2 = list(k); k2[b] -= 1; k2[a] += 1
        if min(k2) < 1 or k2[a] > d - 1: continue
        c2 = boot(adj, k2, 20, 11)
        if c2 >= cur or (c2 >= cur - 1 and hr.random() < 0.3): k, cur = k2, c2; arcs[i] = (b, a)
        if cur > best[0]: best = (cur, list(k))
    k = best[1]; h = {}
    for x in k: h[x] = h.get(x, 0) + 1
    print(f'q={q}: start (1,{q}) Boot {start} -> best visited {boot(adj, k, 200, 3)}; k-hist {sorted(h.items())}', flush=True)
