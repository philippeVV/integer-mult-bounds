import json, collections, itertools
D = json.load(open('geo.json')); n = D['n']; q = D['q']; nt = D['Hnormal']; G = D['G']
dot = lambda a, b: bin(a & b).count('1') & 1
def rank2(vs):
    b = []
    for v in vs:
        for x in b: v = min(v, v ^ x)
        if v: b.append(v)
    return len(b)
c = collections.Counter(); c2 = collections.Counter()
for t, a in itertools.combinations(range(n), 2):
    BK = [s for s in range(n) if not dot(q[s], nt[t]) and not dot(q[s], nt[a])]
    if not BK: c['parallel'] += 1; continue
    plus = [q[s] for s in BK if G[t][s] == G[a][s]]; minus = [q[s] for s in BK if G[t][s] == -G[a][s]]
    assert len(plus) + len(minus) == len(BK)
    c[(len(BK), len(plus), rank2(plus), len(minus), rank2(minus))] += 1
for k, v in sorted(c.items(), key=str): print(k, v)
