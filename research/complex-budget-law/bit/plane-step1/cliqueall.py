"""F2 rank certificate for the plane-masked target: a set Z of ports that is pairwise non-orthogonal (as lines/planes)
forces N[Z,Z] = I_|Z| (diagonal 1, all off-diagonal entries required zero), hence rank_F2 N >= |Z| and #terms >= |Z|.
Z = S plus a maximum pairwise-odd set of lines among C_S = {b live, b not in S : E_b not orthogonal to E_a for all a in S}.
(exact max clique by Bron-Kerbosch with pivoting on <= 111 vertices.)"""
import core, json, sys, time
from collections import Counter
L = list(range(111))
def maxclique(V, adj):
    best = []
    def bk(R, P, X):
        nonlocal best
        if not P and not X:
            if len(R) > len(best): best = R[:]
            return
        if len(R) + len(P) <= len(best): return
        u = max(P | X, key=lambda v: len(P & adj[v]))
        for v in list(P - adj[u]):
            bk(R + [v], P & adj[v], X & adj[v]); P = P - {v}; X = X | {v}
    bk([], set(V), set()); return best
sets = [s for s in json.load(open('sets_all.json')) if len(s[1]) == 5 and s[2] == 4]
# line-only control: max pairwise-odd set among live lines
adjL = {a: {b for b in L if b != a and core.G[a][b] % 2} for a in L}
print('line control omega =', len(maxclique(L, adjL)))
dist = Counter(); t = time.time(); certs = []
for r, S, rk in sets:
    Sset = set(S)
    C = [b for b in L if b not in Sset and all(core.G[b][a] % 2 or core.G[b][r] != 0 for a in S)]
    Q = maxclique(C, {b: adjL[b] & set(C) for b in C})
    Z = list(S) + Q
    # exact check: every pair in Z is non-orthogonal as port spaces
    E = {a: [core.LINES[a]] + ([core.LINES[r]] if a in Sset else []) for a in Z}
    assert all(any(core.dot(u, v) for u in E[a] for v in E[b]) for a in Z for b in Z if a != b)
    dist[len(Z)] += 1; certs.append([r, S, Q])
print('plane omega lower-bound certificates (|Z| distribution):', dict(dist), round(time.time() - t), 's')
json.dump(certs, open('clique_certs.json', 'w'))
