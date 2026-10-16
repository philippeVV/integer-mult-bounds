"""Bootstrap order search with tie-break 'place deficiency-0 vertices as soon as j = k' (cost = sum (j-k)^+ too),
returns best order; then checks the back-arc flow (schedule from the order: Eulerian orientation whose
back arcs number exactly cost(order)).  usage: boot3.py GRAPH restarts [expected_cost] (graph file via env GF)"""
import sys, json, random, os
def adjacency(G):
    V = G['Pts'] + G['Lns']; E = G['E']; ends = [[] for _ in range(E)]
    for v, es in enumerate(V):
        for e in es: ends[e].append(v)
    adj = [[] for _ in V]
    for a, b in ends: adj[a].append(b); adj[b].append(a)
    return adj, [len(es) // 2 for es in V], ends
def greedy(adj, k, rng):
    n = len(adj); j = [0] * n; placed = [False] * n; order = []; c = 0
    # always place a vertex with j >= k if one exists (prefer the smallest j), else max j
    import heapq
    ready = []  # (j, tiebreak, v) for j >= k
    cnt = {}
    pool = set(range(n))
    byj = [set() for _ in range(max(len(a) for a in adj) + 1)]
    for v in range(n): byj[0].add(v)
    for _ in range(n):
        v = None
        for jj in range(len(byj)):
            if jj >= min(k) and byj[jj]:
                cands = [w for w in byj[jj] if jj >= k[w]]
                if cands: v = rng.choice(cands); break
        if v is None:
            for jj in range(len(byj) - 1, -1, -1):
                if byj[jj]: v = rng.choice(list(byj[jj])); break
        byj[j[v]].discard(v); placed[v] = True; order.append(v); c += max(k[v] - j[v], 0)
        for w in adj[v]:
            if not placed[w]:
                byj[j[w]].discard(w); j[w] += 1; byj[j[w]].add(w)
    return c, order
def flow_ok(adj, k, order):
    """excess vertices (j>k) must send j-k back arcs to earlier nbrs with deficit (k-j>0); max-flow check."""
    pos = {v: i for i, v in enumerate(order)}; n = len(adj)
    jv = [sum(1 for w in adj[v] if pos[w] < pos[v]) for v in range(n)]
    import networkx as nx
    Gf = nx.DiGraph(); src, snk = 's', 't'; need = 0
    for v in range(n):
        if jv[v] > k[v]: Gf.add_edge(src, ('o', v), capacity=jv[v] - k[v]); need += jv[v] - k[v]
        if jv[v] < k[v]: Gf.add_edge(('i', v), snk, capacity=k[v] - jv[v])
        for w in adj[v]:
            if pos[w] < pos[v] and jv[v] > k[v] and jv[w] < k[w]:
                Gf.add_edge(('o', v), ('i', w), capacity=1)
    if need == 0: return True, 0, 0
    val = nx.maximum_flow_value(Gf, src, snk)
    return val == need, val, need
if __name__ == '__main__':
    G = json.load(open(os.environ.get('GF', 'graphs.json')))[sys.argv[1]]; R = int(sys.argv[2])
    adj, k, ends = adjacency(G); rng = random.Random(3); best = None
    for _ in range(R):
        c, o = greedy(adj, k, rng)
        if best is None or c < best[0]: best = (c, o)
    ok = flow_ok(adj, k, best[1])
    print(json.dumps(dict(graph=sys.argv[1], E=G['E'], boot=best[0], s=round(best[0] / G['E'], 4), schedule_flow=ok)))
    if len(sys.argv) > 3:   # limits-1 replay: expected greedy cost (the flow verdict is reported, not asserted)
        assert best[0] == int(sys.argv[3]), best[0]
        print('PASS_EXPXE_BOOT3')
