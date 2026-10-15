import json, collections, itertools
D = json.load(open('geo.json')); n = D['n']; q = D['q']; nt = D['Hnormal']
dot = lambda a, b: bin(a & b).count('1') & 1
src_of = {q[s]: s for s in range(n)}
print('q==normal?', sum(q[t] == nt[t] for t in range(n)), ' normals that are labels', sum(nt[t] in src_of for t in range(n)))
nset = {nt[t]: t for t in range(n)}
trip = sum(1 for a, b in itertools.combinations(range(n), 2) if (nt[a] ^ nt[b]) in nset) // 1
print('pairs (a,b) whose third hyperplane is also a target', trip, ' (each triple counted 3x) triples', trip // 3)
inK = lambda v, t, w: (not dot(v, nt[t])) and (not dot(v, w))
# sources per K = H_a ∩ H_b
cnt = collections.Counter()
for a, b in itertools.combinations(range(n), 2):
    cnt[sum(1 for s in range(n) if not dot(q[s], nt[a]) and not dot(q[s], nt[b]))] += 1
print('sources in H_a∩H_b distribution', sorted(cnt.items()))
# for each target: hyperplanes K of H_t (given by w mod n_t) ; valid cover pairs
ncov = []; ncov_tt = []
for t in range(n):
    S_ = [q[s] for s in range(n) if not dot(q[s], nt[t])]
    ws = [w for w in range(1, 512) if w != nt[t] and w < (w ^ nt[t])]
    good = 0; goodtt = 0
    for w1, w2 in itertools.combinations(ws, 2):
        if all(not dot(v, w1) or not dot(v, w2) for v in S_):
            good += 1
            if ((w1 in nset or (w1 ^ nt[t]) in nset) and (w2 in nset or (w2 ^ nt[t]) in nset)): goodtt += 1
    ncov.append(good); ncov_tt.append(goodtt)
print('valid 2-covers per target', collections.Counter(ncov))
print('valid 2-covers with both K = H_t ∩ H_other', collections.Counter(ncov_tt))
