"""Fast bootstrap bound with bucket queue; supports per-side k (kP, kL) for non-half-rate checks.
usage: boot2.py GRAPH restarts [kP kL]   (default k = deg/2)"""
import sys, json, random
def run(adj, k, rng):
    n = len(adj); j = [0] * n; placed = [False] * n
    K = max(k); buckets = [[] for _ in range(K + 1)]; where = [None] * n
    def put(w):
        d = max(k[w] - j[w], 0); b = buckets[d]; where[w] = (d, len(b)); b.append(w)
    def take(w):
        d, i = where[w]; b = buckets[d]; last = b[-1]; b[i] = last; where[last] = (d, i); b.pop()
    for w in range(n): put(w)
    c = 0
    for _ in range(n):
        d = next(d for d in range(K + 1) if buckets[d])
        b = buckets[d]; v = b[rng.randrange(len(b))]; take(v); placed[v] = True; c += d
        for w in adj[v]:
            if not placed[w]:
                take(w); j[w] += 1; put(w)
    return c
if __name__ == '__main__':
    import os; G = json.load(open(os.environ.get('GF', 'graphs.json')))[sys.argv[1]]; R = int(sys.argv[2])
    V = G['Pts'] + G['Lns']; E = G['E']; nP = len(G['Pts'])
    ends = [[] for _ in range(E)]
    for v, es in enumerate(V):
        for e in es: ends[e].append(v)
    adj = [[] for _ in V]
    for a, b in ends: adj[a].append(b); adj[b].append(a)
    if len(sys.argv) > 4: k = [int(sys.argv[3])] * nP + [int(sys.argv[4])] * (len(V) - nP)
    else: k = [len(es) // 2 for es in V]
    rng = random.Random(7); best = min(run(adj, k, rng) for _ in range(R))
    print(json.dumps(dict(graph=sys.argv[1], E=E, girth=G['girth'], k=[k[0], k[-1]], boot=best, s_bound=round(best / E, 4))))
