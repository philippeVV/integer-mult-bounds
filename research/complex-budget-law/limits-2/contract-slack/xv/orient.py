"""Best schedulable Eulerian orientation for a graph from expXP graphs.json.
Order sigma: best of (expXP bkey, greedy restarts + annealing via expXP/bootsa).  For that order, the Eulerian (k-in/k-out)
orientation with the fewest back arcs F(sigma) by min-cost flow (networkx).  R = back arcs; Gamma - R is acyclic (forward in sigma).
usage: orient.py NAME restarts sweeps  -> orient_NAME.json  {order, tail[E], head[E], R, F, cost}"""
import sys, json, os, numpy as np, networkx as nx
HERE = os.path.dirname(os.path.abspath(__file__)); XP = os.path.join(HERE, '..', 'expXP'); sys.path.insert(0, XP)
import bootsa as bs

name, R, sweeps = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
G, n, E, ptr, adj, deg = bs.load(os.path.join(XP, 'graphs.json'), name)
k = deg // 2
EV = G['Pts'] + G['Lns']
ends = [[] for _ in range(E)]
for v, es in enumerate(EV):
    for e in es: ends[e].append(v)
ends = np.array(ends)

def F_of(key):
    pos = np.argsort(np.argsort(key))
    a, b = ends[:, 0], ends[:, 1]; fw = pos[a] < pos[b]
    u = np.where(fw, a, b); w = np.where(fw, b, a)          # forward arc u -> w
    j = np.bincount(w, minlength=n); cost = int(np.maximum(k - j, 0).sum())
    D = nx.DiGraph()
    for v in range(n): D.add_node(v, demand=int(k[v] - j[v]))   # networkx: demand = inflow - outflow
    # flipping edge e (w->u) moves one unit of in-degree from w to u: flow w -> u
    M = nx.MultiDiGraph(D)
    for e in range(E): M.add_edge(int(w[e]), int(u[e]), key=e, capacity=1, weight=1)
    fc, fd = nx.network_simplex(M)
    flip = np.zeros(E, bool)
    for x in fd:
        for y in fd[x]:
            for e, val in fd[x][y].items():
                if val: flip[e] = True
    tail = np.where(flip, w, u); head = np.where(flip, u, w)
    back = flip.copy()
    indeg = np.bincount(head, minlength=n); assert (indeg == k).all()
    return cost, int(fc), tail, head, np.flatnonzero(back), pos

rng = np.random.default_rng(0)
cands = [np.load(os.path.join(XP, f'bkey_{name}.npy'))]
best = None
for r in range(R):
    key = bs.greedy_key(n, ptr, adj, k, rng)
    c, key2 = bs.anneal(key, ptr, adj, k, sweeps, 0.6, 0.05, 1000 + r)
    cands.append(bs.renorm(key2))
for key in cands:
    c, Fv, tail, head, back, pos = F_of(key)
    if best is None or (Fv, c) < (best[1], best[0]): best = (c, Fv, tail, head, back, pos)
c, Fv, tail, head, back, pos = best
# acyclicity check of Gamma - R (Kahn)
Dg = nx.DiGraph(); Dg.add_nodes_from(range(n))
for e in range(E):
    if e not in set(back.tolist()): Dg.add_edge(int(tail[e]), int(head[e]))
assert nx.is_directed_acyclic_graph(Dg)
out = dict(name=name, n=n, E=E, cost=c, F=Fv, order=[int(x) for x in np.argsort(pos)], tail=tail.tolist(), head=head.tolist(),
           R=back.tolist(), E_over_6=E / 6)
json.dump(out, open(os.path.join(HERE, f'orient_{name}.json'), 'w'))
print(json.dumps({x: out[x] for x in ['name', 'n', 'E', 'cost', 'F', 'E_over_6']}))
