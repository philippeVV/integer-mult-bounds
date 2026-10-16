"""PG(2,16) full, profile (8 pts, 9 lines): best Boot order (anneal) and the exact minimum number of back
arcs B(sigma) over all orientations with in-degrees k (min-cost flow).  B(sigma) = |R| of an actual causal
schedule (R = back arcs, Gamma-R acyclic by sigma).  |R| = dim + c needs c zero-source and c erased edges
(free under idea 2); the a=1/9 threshold needs (900-c)/(4641-2c) > 1/(1+5^(8/9)), i.e. c <= 6.
usage: sched16.py iters seed   (needs networkx)"""
import sys, os, random, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault('T0', '0.7')
import networkx as nx
from xw import pg, greedy, cost
from boot_anneal import anneal

def min_back_arcs(adj, edges, k, order):
    pos = {v: i for i, v in enumerate(order)}
    G = nx.DiGraph(); E = len(edges)
    G.add_node('s', demand=-E)
    for v in range(len(adj)): G.add_edge(('v', v), 't', capacity=k[v], weight=0)
    G.add_node('t', demand=E)
    for i, (a, b) in enumerate(edges):
        u, w = (a, b) if pos[a] < pos[b] else (b, a)      # u earlier
        G.add_edge('s', ('e', i), capacity=1, weight=0)
        G.add_edge(('e', i), ('v', w), capacity=1, weight=0)   # forward arc u->w
        G.add_edge(('e', i), ('v', u), capacity=1, weight=1)   # back arc w->u
    f = nx.min_cost_flow(G)
    return sum(f[('e', i)][('v', min(e, key=lambda x: pos[x]))] for i, e in enumerate(edges))

if __name__ == '__main__':
    iters, seed = int(sys.argv[1]), int(sys.argv[2])
    N, edges, adj = pg(16); E = len(edges); k = [8] * N + [9] * N
    rng = random.Random(seed)
    g = min((greedy(adj, k, rng) for _ in range(200)), key=lambda t: t[0])
    c, o = anneal(adj, k, iters, rng, g[1])
    if g[0] <= c: c, o = g
    B = min_back_arcs(adj, edges, k, o)
    thr = 1 / (1 + 5 ** (8 / 9)); cc = B - 900
    print(json.dumps(dict(boot_upper=c, min_back_arcs_this_order=B, c=cc,
                          s_after_CJ=round((900 - cc) / (E - 2 * cc), 5), threshold=round(thr, 5))))
