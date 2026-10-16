"""Annealed bootstrap-order search (upper bounds on Boot(Gamma,k)) + back-arc transport check.

Order = real keys; a move re-keys one vertex (near a random neighbour, or uniformly); the cost
sum_v (k_v - j_v)^+ changes only at v and its neighbours, so a move is O(deg).
usage: boot_anneal.py q kP kL iters seed      (PG(2,q) Singer incidence from xw.py)
Prints best cost, and whether some orientation with in-degrees k has exactly that many back arcs
(the schedule condition of expXE section 2.3, generalised to any in-degree profile)."""
import sys, random, math, json
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from xw import pg, greedy, cost

def transport(adj, k, order):
    """Max-flow: excess vertices (j>k) send j-k back arcs to earlier deficit nbrs (k>j). Simple
    augmenting paths on a unit bipartite network."""
    pos = [0] * len(adj)
    for i, v in enumerate(order): pos[v] = i
    j = [sum(1 for w in adj[v] if pos[w] < pos[v]) for v in range(len(adj))]
    supply = {v: j[v] - k[v] for v in range(len(adj)) if j[v] > k[v]}
    demand = {v: k[v] - j[v] for v in range(len(adj)) if j[v] < k[v]}
    need = sum(supply.values()); assert need == sum(demand.values())
    # bipartite b-matching on edges (v excess, u deficit, u earlier nbr of v); augmenting paths
    used = {}                       # (v,u) -> 1 if chosen
    load_u = {u: 0 for u in demand}; load_v = {v: 0 for v in supply}
    nbr = {v: [u for u in adj[v] if u in demand and pos[u] < pos[v]] for v in supply}
    def augment(v, seen):
        for u in nbr[v]:
            if used.get((v, u)) or u in seen: continue
            seen.add(u)
            if load_u[u] < demand[u]:
                used[(v, u)] = 1; load_u[u] += 1; return True
            for v2 in [x for (x, y), f in used.items() if y == u and f]:
                if augment(v2, seen):
                    used[(v2, u)] = 0; used[(v, u)] = 1; return True
        return False
    flow = 0
    for v in supply:
        for _ in range(supply[v]):
            if augment(v, set()): flow += 1
    return flow, need

def anneal(adj, k, iters, rng, start):
    n = len(adj); key = [0.0] * n
    for i, v in enumerate(start): key[v] = (i + rng.random()) / n
    j = [sum(1 for w in adj[v] if key[w] < key[v]) for v in range(n)]
    c = sum(max(k[v] - j[v], 0) for v in range(n)); best = c; bestkey = list(key)
    T0 = float(__import__("os").environ.get("T0", "2.0"))
    for it in range(iters):
        T = T0 * (1 - it / iters) + 0.02
        v = rng.randrange(n)
        if rng.random() < 0.8:
            w = rng.choice(adj[v]); nk = key[w] + (rng.random() - 0.5) * 2e-4
        else:
            nk = rng.random()
        old = key[v]; delta = 0; dj = 0; changes = []
        for w in adj[v]:
            before_old = key[w] < old; before_new = key[w] < nk
            if before_old != before_new:
                # w's relation flips: if w was before v (j_v counted it) now after -> j_v -1, j_w +1
                s = -1 if before_old else 1
                dj += s
                jw = j[w]; nw = jw - s
                delta += max(k[w] - nw, 0) - max(k[w] - jw, 0); changes.append((w, nw))
        nv = j[v] + dj; delta += max(k[v] - nv, 0) - max(k[v] - j[v], 0)
        if delta <= 0 or rng.random() < math.exp(-delta / T):
            key[v] = nk; j[v] = nv
            for w, nw in changes: j[w] = nw
            c += delta
            if c < best: best = c; bestkey = list(key)
    order = sorted(range(n), key=lambda x: bestkey[x])
    assert cost(adj, k, order) == best
    return best, order

if __name__ == '__main__':
    q, kP, kL, iters, seed = map(int, sys.argv[1:6])
    N, edges, adj = pg(q); E = len(edges); k = [kP] * N + [kL] * N
    assert sum(k) == E
    rng = random.Random(seed)
    g = min((greedy(adj, k, rng) for _ in range(200)), key=lambda t: t[0])
    print(json.dumps(dict(q=q, E=E, greedy=g[0])), flush=True)
    best = g
    for r in range(3):
        c, o = anneal(adj, k, iters, rng, best[1])
        if c < best[0]: best = (c, o)
        print(json.dumps(dict(round=r, anneal=c, best=best[0], s=round(best[0] / E, 5))), flush=True)
    fl, need = transport(adj, k, best[1])
    print(json.dumps(dict(q=q, k=[kP, kL], E=E, boot_upper=best[0], s_upper=round(best[0] / E, 5),
                          transport_flow=fl, transport_need=need, schedulable_order=fl == need)))
