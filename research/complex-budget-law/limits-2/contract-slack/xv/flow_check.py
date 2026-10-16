"""expXV flow-screen check on PG(2,5) (numpy, scipy/HiGHS, networkx; a few seconds).
1. orient_PG5.json is an Eulerian 3-in/3-out orientation of PG(2,5) with |R| = F = 37 back arcs and Gamma - R acyclic
   (Kahn), so for in/out-invertible gates dim W <= |R| = 37.
2. Flow screen at S = 32 sessions drawn from R (as flowtest.py): every session has distance >= 5 in Gamma - Rs (girth 6), so
   F_* <= (E - S)/5 = 30.8 < S; the exact undirected multiflow LP attains F_* = 30.8, F_*/W = 0.9625.  So F_* < W is
   automatic once S > E/6 and the flow screen is never decisive on these hosts (README, "Screens").
Prints PASS_EXPXV_FLOW."""
import os, sys, json
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import numpy as np
from flowlp import fstar
G = json.load(open(os.path.join(H, '..', '..', '..', 'limits-1', 'complex-dual', 'xg', 'graphs_small.json')))['PG5']
EV = G['Pts'] + G['Lns']; E = G['E']; n = len(EV)
ends = [[] for _ in range(E)]
for v, es in enumerate(EV):
    for e in es: ends[e].append(v)
ends = np.array(ends); OR = json.load(open(os.path.join(H, 'orient_PG5.json')))
tail, head, R = OR['tail'], OR['head'], set(OR['R'])
ok = all(sorted((tail[e], head[e])) == sorted(ends[e].tolist()) for e in range(E))   # arc e orients edge e
indeg = np.bincount(head, minlength=n); outdeg = np.bincount(tail, minlength=n)
ok &= set(indeg) == {3} and set(outdeg) == {3} and len(R) == OR['F'] == 37
adj = [[] for _ in range(n)]; deg = [0] * n
for e in range(E):
    if e not in R: adj[tail[e]].append(head[e]); deg[head[e]] += 1
st = [v for v in range(n) if deg[v] == 0]; seen = 0
while st:
    v = st.pop(); seen += 1
    for w in adj[v]:
        deg[w] -= 1
        if deg[w] == 0: st.append(w)
ok &= seen == n
print(f'orientation: Eulerian 3/3 {set(indeg) == {3}}, |R| = {len(R)} (F = {OR["F"]}), Gamma - R acyclic {seen == n}')
rng = np.random.default_rng(0)
for t in range(2):
    res = fstar(ends, rng.choice(sorted(R), 32, replace=False), n); print(json.dumps(res))
    ok &= res['status'] == 0 and res['dist_min'] >= 5 and abs(res['F_star'] - res['volume_bound']) < 1e-6 and res['F_over_W'] < 0.98
print('PASS_EXPXV_FLOW' if ok else 'FAIL_EXPXV_FLOW')
