import json, collections, itertools
from ortools.sat.python import cp_model
D = json.load(open('geo.json')); n = D['n']; q = D['q']; nt = D['Hnormal']
dot = lambda a, b: bin(a & b).count('1') & 1
def span(bs):
    s = {0}
    for b in bs:
        if b not in s: s |= {x ^ b for x in s}
    return frozenset(s)
Hs = [frozenset(v for v in range(512) if not dot(v, nt[t])) for t in range(n)]
supp = {t: [s for s in range(n) if q[s] in Hs[t]] for t in range(n)}
assert all(len(v) == 48 for v in supp.values())
# later-frame (non-line) access: does X_s have a retained frame of dim>=2 inside H_t?
later = {(t, s): any(seg['dim'] >= 2 and span(seg['basis']) <= Hs[t] for seg in D['X'][s]) for t in range(n) for s in supp[t]}
print('pairs with later frame inside', sum(later.values()))
Ts = {s: [t for t in range(n) if q[s] in Hs[t]] for s in range(n)}
Tline = {s: [t for t in Ts[s] if not later[t, s]] for s in range(n)}
def mincover(univ, sets):
    m = cp_model.CpModel(); x = [m.NewBoolVar('') for _ in sets]
    for u in univ: m.AddBoolOr([x[i] for i, S in enumerate(sets) if u in S])
    m.Minimize(sum(x)); sol = cp_model.CpSolver(); sol.parameters.num_workers = 1; sol.parameters.max_time_in_seconds = 60
    r = sol.Solve(m); assert r == cp_model.OPTIMAL, r
    return int(sol.ObjectiveValue()), [i for i in range(len(sets)) if sol.Value(x[i])]
kap = {}; kap_all = {}
for s in range(n):
    ws = [w for w in range(1, 512) if w != q[s] and w < (w ^ q[s])]
    sets = [frozenset(t for t in Ts[s] if not dot(nt[t], w)) for w in ws]
    kap[s] = mincover(Tline[s], sets)[0]
    kap_all[s] = mincover(Ts[s], sets)[0]
print('kappa_s (line-needed targets) distribution', sorted(collections.Counter(kap.values()).items()), 'sum', sum(kap.values()))
print('kappa_s (all 48 targets) distribution', sorted(collections.Counter(kap_all.values()).items()), 'sum', sum(kap_all.values()))
print('Tline sizes', sorted(collections.Counter(len(v) for v in Tline.values()).items()))
# rank of normals within q_s-perp
def rank2(vs):
    b = []
    for v in vs:
        for x in b: v = min(v, v ^ x)
        if v: b.append(v)
    return len(b)
print('rank of needed normals per source', collections.Counter(rank2([nt[t] for t in Ts[s]]) for s in range(n)))
