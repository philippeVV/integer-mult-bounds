"""Bootstrap ratio cost/E for random d-regular bipartite graphs (union of d random matchings), n per side."""
import sys, random, json
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from boot2 import run
d, n, R = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
rng = random.Random(d * 1000 + n)
adj = [[] for _ in range(2 * n)]; seen = set()
for _ in range(d):
    while True:
        p = list(range(n)); rng.shuffle(p)
        if all((i, p[i]) not in seen for i in range(n)): break
    for i in range(n):
        seen.add((i, p[i])); adj[i].append(n + p[i]); adj[n + p[i]].append(i)
k = [d // 2] * (2 * n); E = d * n
best = min(run(adj, k, rng) for _ in range(R))
print(json.dumps(dict(d=d, n=n, E=E, boot=best, ratio=round(best / E, 4))))
print('PASS_EXPXE_BOOTRAND')
