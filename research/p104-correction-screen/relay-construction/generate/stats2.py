import json, collections, itertools
D = json.load(open('geo.json')); n = D['n']; q = D['q']; nt = D['Hnormal']
dot = lambda a, b: bin(a & b).count('1') & 1
def sp(bs):
    s = {0}
    for b in bs:
        if b not in s: s |= {x ^ b for x in s}
    return frozenset(s)
def perp(ws): return frozenset(v for v in range(512) if all(not dot(v, w) for w in ws))
def rank2(vs):
    b = []
    for v in vs:
        for x in b: v = min(v, v ^ x)
        if v: b.append(v)
    return len(b)
supp = {t: [s for s in range(n) if not dot(q[s], nt[t])] for t in range(n)}
U = {}
for t in range(n):
    us = [u for u in range(1, 512) if all(dot(u, q[s]) for s in supp[t])]
    U[t] = us
print('u_t candidates per target', collections.Counter(len(v) for v in U.values()))
print('u_t == q_t possible', sum(q[t] in U[t] for t in range(n)))
# canonical u_t mod n_t: two candidates u, u^n_t
u = {t: min(U[t]) for t in range(n)}
print('u candidates are {u,u^n}:', all(set(U[t]) == {u[t], u[t] ^ nt[t]} for t in range(n)))
# Does the binary pairing <q_t, q_s> separate supp?
print('pairing <q_t,q_s>=1 on supp, count of t', sum(all(dot(q[t], q[s]) for s in supp[t]) for t in range(n)))
# J0 for pairs
J0s = collections.defaultdict(list)
for t, a in itertools.combinations(range(n), 2):
    J = perp([nt[t], nt[a], u[t] ^ u[a]])
    BK = [s for s in range(n) if not dot(q[s], nt[t]) and not dot(q[s], nt[a])]
    assert all(q[s] in J for s in BK)
    J0s[J].append((t, a, len(BK)))
print('distinct J0', len(J0s), 'dims', collections.Counter(rank2(sorted(J)) for J in J0s) )
print('pairs per J0', collections.Counter(len(v) for v in J0s.values()))
