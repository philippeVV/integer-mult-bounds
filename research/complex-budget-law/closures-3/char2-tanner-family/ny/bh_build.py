"""Exact (stdlib) reconstruction of the Braverman-He PG(2,9) causal code (arXiv:2610.10108 Sec. 4),
plus the reversible-gate reading needed by the repo's bit contract (04-swap.tex:105-123).

Checks (all exact, F9 = F3[t]/(t^2+1), element a+b*t encoded as a+3b, matching eta_r = r):
  1. 182 vertices, 910 candidate edges, 10-regular, no 4-cycles; 735 physical edges, 157 sessions.
  2. orientation: d-(v)=d+(v)=5, d_C^+(v)=d_J^-(v) at every vertex.
  3. each vertex's local [10,5] code is MDS (every 5 positions invertible)  -> the vertex is a
     reversible 5-in/5-out F9-linear gate K(v) -> U(v).
  4. causal simulation: all 157 sinks decode their source symbol for 157 unit inputs + random inputs.
  5. dist_G(s_q,t_q) >= 5 for all q and sum of distances (routing bound 735/sum).
  6. reversible reading: 18 C-ports become arbitrary helper inputs, 18 J-ends become outputs;
     global map T on F9^175 = [[I + ?]]; we report the blocks L (sink contamination by helpers),
     G, H and ranks, i.e. what a permutation completion has to undo.
Writes bh_graph.json for the flow LPs.
"""
import json, random, sys
from collections import deque
from bh_data import R, PI, C, J

# ---------- F9 ----------
def f9(a, b): return (a % 3) + 3 * (b % 3)
def ab(x): return x % 3, x // 3
ADD = [[f9(ab(x)[0] + ab(y)[0], ab(x)[1] + ab(y)[1]) for y in range(9)] for x in range(9)]
MUL = [[f9(ab(x)[0]*ab(y)[0] - ab(x)[1]*ab(y)[1], ab(x)[0]*ab(y)[1] + ab(x)[1]*ab(y)[0]) for y in range(9)] for x in range(9)]
NEG = [f9(-ab(x)[0], -ab(x)[1]) for x in range(9)]
INV = [None] + [next(y for y in range(9) if MUL[x][y] == 1) for x in range(1, 9)]
def sub(x, y): return ADD[x][NEG[y]]
def pw(x, k):
    r = 1
    for _ in range(k): r = MUL[r][x]
    return r
assert all(MUL[x][INV[x]] == 1 for x in range(1, 9))
assert MUL[3][3] == NEG[1]  # t^2 = -1

def solve(A, b):
    """Solve square A x = b over F9 (A invertible) by Gauss-Jordan; returns None if singular."""
    n = len(A); M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c]), None)
        if p is None: return None
        M[c], M[p] = M[p], M[c]
        iv = INV[M[c][c]]; M[c] = [MUL[iv][v] for v in M[c]]
        for r in range(n):
            if r != c and M[r][c]:
                f = M[r][c]; M[r] = [sub(M[r][k], MUL[f][M[c][k]]) for k in range(n + 1)]
    return [M[i][n] for i in range(n)]

def rank(A):
    M = [r[:] for r in A]; rk = 0; ncol = len(M[0]) if M else 0
    for c in range(ncol):
        p = next((r for r in range(rk, len(M)) if M[r][c]), None)
        if p is None: continue
        M[rk], M[p] = M[p], M[rk]
        iv = INV[M[rk][c]]; M[rk] = [MUL[iv][v] for v in M[rk]]
        for r in range(len(M)):
            if r != rk and M[r][c]:
                f = M[r][c]; M[r] = [sub(M[r][k], MUL[f][M[rk][k]]) for k in range(ncol)]
        rk += 1
    return rk

# ---------- PG(2,9) ----------
T = [(0, 0, 1)] + [(0, 1, r) for r in range(9)] + [(1, r, s) for r in range(9) for s in range(9)]
assert len(T) == 91
def dot(p, l): return ADD[ADD[MUL[p[0]][l[0]]][MUL[p[1]][l[1]]]][MUL[p[2]][l[2]]]
nbrs = {i: sorted(j for j in range(91) if dot(T[i], T[j]) == 0) for i in range(91)}
assert all(len(v) == 10 for v in nbrs.values())
edges = []  # e_q = (point i, line vertex 91+j)
for i in range(91):
    for j in nbrs[i]: edges.append((i, 91 + j))
assert len(edges) == 910
inc = {v: [] for v in range(182)}
for q, (a, b) in enumerate(edges): inc[a].append(q); inc[b].append(q)
assert all(len(inc[v]) == 10 for v in range(182))
# no 4-cycles: two points share <=1 line
for i in range(91):
    for k in range(i + 1, 91):
        assert len(set(nbrs[i]) & set(nbrs[k])) == 1

Rs, Cs, Js = set(R), set(C), set(J)
assert Cs <= Rs and Js <= Rs and not (Cs & Js)
B = sorted(Rs - Cs - Js); assert len(B) == 157
phys = [q for q in range(910) if q not in Rs]; assert len(phys) == 735
pos = {v: k for k, v in enumerate(PI)}

# orientation
def tail_head(q):
    a, b = edges[q]
    early, late = (a, b) if pos[a] < pos[b] else (b, a)
    return (early, late) if q not in Rs else (late, early)
din = {v: [] for v in range(182)}; dout = {v: [] for v in range(182)}
for q in range(910):
    t, h = tail_head(q); dout[t].append(q); din[h].append(q)
for v in range(182):
    assert len(din[v]) == 5 and len(dout[v]) == 5, v
    assert sum(q in Cs for q in dout[v]) == sum(q in Js for q in din[v])

# local codes: generator rows (evaluation vectors) per incident edge
def coords(z, y):
    perp = [t for t in T if dot(z, t) == 0]
    u, w = perp[0], perp[1]
    for (c1, c2) in [(0, 1), (0, 2), (1, 2)]:
        sol = solve([[u[c1], w[c1]], [u[c2], w[c2]]], [y[c1], y[c2]])
        if sol is not None:
            lam, mu = sol
            assert all(ADD[MUL[lam][u[k]]][MUL[mu][w[k]]] == y[k] for k in range(3))
            return lam, mu
    raise AssertionError
gen = {}
for v in range(182):
    z = T[v] if v < 91 else T[v - 91]
    for q in inc[v]:
        other = edges[q][1] - 91 if v < 91 else edges[q][0]
        lam, mu = coords(z, T[other])
        if v >= 91: lam, mu = pw(lam, 3), pw(mu, 3)
        gen[(v, q)] = [MUL[pw(lam, 4 - r)][pw(mu, r)] for r in range(5)]
# MDS check: every 5-subset of the 10 rows invertible (sufficient: check all 252 subsets per vertex)
from itertools import combinations
for v in range(182):
    rows = [gen[(v, q)] for q in inc[v]]
    for S in combinations(range(10), 5):
        assert rank([rows[k] for k in S]) == 5
print("local codes: all 182 vertices MDS [10,5] over F9 (252 subsets each)")

sess_src = {q: tail_head(q)[1] for q in B}   # B-edge oriented later->earlier: head = earlier = s_q
sess_snk = {q: tail_head(q)[0] for q in B}
for q in B:
    assert pos[sess_src[q]] < pos[sess_snk[q]]

def run(x, cin):
    """x: dict q->F9 source symbol (q in B); cin: dict (q, v)->value for C-edge ends (0 = BH protocol).
    Returns (sink outputs dict, J-end outputs dict, all edge values seen)."""
    val = {}          # physical edge values
    out_sink, out_j = {}, {}
    for v in PI:
        known = []
        for q in din[v]:
            if q in Cs: known.append((q, cin.get((q, v), 0)))
            elif q in Js: pass
            elif q in Rs: known.append((q, x[q]))       # incoming B edge: own source
            else: known.append((q, val[q]))
        for q in dout[v]:
            if q in Cs: known.append((q, cin.get((q, v), 0)))
        assert len(known) == 5
        coef = solve([gen[(v, q)] for q, _ in known], [c for _, c in known])
        def ev(q):
            g = gen[(v, q)]; s = 0
            for r in range(5): s = ADD[s][MUL[g[r]][coef[r]]]
            return s
        for q in dout[v]:
            if q in Cs: continue
            if q in Js: out_j[(q, v)] = ev(q)
            elif q in Rs: out_sink[q] = ev(q)
            else: val[q] = ev(q)
        for q in din[v]:
            if q in Js: out_j[(q, v)] = ev(q)
    return out_sink, out_j

# 4. correctness of the causal code (linear => unit vectors + random suffice; we do both)
zero = {q: 0 for q in B}
for q in B:
    x = dict(zero); x[q] = 1
    o, _ = run(x, {})
    assert o == x, q
rng = random.Random(1)
for _ in range(5):
    x = {q: rng.randrange(9) for q in B}
    assert run(x, {})[0] == x
print("causal code: 157 sinks decode correctly (157 unit + 5 random inputs)")

# 5. distances in physical graph
adj = {v: [] for v in range(182)}
for q in phys:
    a, b = edges[q]; adj[a].append(b); adj[b].append(a)
def bfs(s):
    d = {s: 0}; dq = deque([s])
    while dq:
        u = dq.popleft()
        for w in adj[u]:
            if w not in d: d[w] = d[u] + 1; dq.append(w)
    return d
dist = {q: bfs(sess_src[q]).get(sess_snk[q]) for q in B}
from collections import Counter
print("session distance histogram:", sorted(Counter(dist.values()).items()), "sum", sum(dist.values()))
print("routing bound 735/sum =", 735, "/", sum(dist.values()))

# 6. reversible reading: helper lanes
cports = [(q, v) for v in range(182) for q in din[v] + dout[v] if q in Cs]
jports = [(q, v) for v in range(182) for q in din[v] + dout[v] if q in Js]
assert len(cports) == 18 and len(jports) == 18
inputs = [('x', q) for q in B] + [('c', p) for p in cports]
def column(lane):
    x = dict(zero); cin = {}
    if lane[0] == 'x': x[lane[1]] = 1
    else: cin[lane[1]] = 1
    o, j = run(x, cin)
    return [o[q] for q in B] + [j[p] for p in jports]
cols = [column(l) for l in inputs]
Tm = [[cols[c][r] for c in range(175)] for r in range(175)]
assert rank(Tm) == 175
L = [row[157:] for row in Tm[:157]]; G = [row[:157] for row in Tm[157:]]; H = [row[157:] for row in Tm[157:]]
assert all(Tm[r][c] == (1 if r == c else 0) for r in range(157) for c in range(157))
nzL = sum(1 for row in L if any(row))
print("reversible map T on F9^175: x-block of sinks = I;  rank L =", rank(L), "(sinks contaminated:", nzL, ")",
      " rank G =", rank(G), " rank H =", rank(H))
json.dump({"edges": [list(edges[q]) for q in phys],
           "sessions": [[sess_src[q], sess_snk[q]] for q in B],
           "cports": [[q, v] for q, v in cports], "jports": [[q, v] for q, v in jports],
           "dist": [dist[q] for q in B],
           "L_rows_nonzero": nzL, "T": Tm,
           "src_vertex": [sess_src[q] for q in B], "snk_vertex": [sess_snk[q] for q in B]}, open("bh_graph.json", "w"))
print("supp sizes: G", sum(1 for r in G for v in r if v), " L", sum(1 for r in L for v in r if v), " H", sum(1 for r in H for v in r if v))
print("wrote bh_graph.json")
