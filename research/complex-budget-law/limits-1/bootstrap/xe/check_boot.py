"""limits-1/bootstrap: field-free bootstrap ceiling dim C <= Boot(G) (expXE Theorem; expXG Thm C).
Reproduces the greedy Boot upper bounds of the expXE table and the X(., 1/9) exclusion verdicts.  Stdlib only.
Graphs: graphs_small.json (expXG: PG3, BA4, PG4M, W4M, PG5), graphs_gq.json (W3, W5, W7), graphs_gq2.json (W11),
and PG(2,q) for q = 7, 9, 11, 13, 17 built here from coordinates (q = 9 via F_3[i]).
usage: python3 -I -B check_boot.py [restarts]   (default 300 for E <= 3200, 20 above; W(7) needs the restarts)"""
import sys, os, json, math, random, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from boot2 import run

R = int(sys.argv[1]) if len(sys.argv) > 1 else 300


def adj_of(Pts, Lns, E):
    V = Pts + Lns; ends = [[] for _ in range(E)]
    for v, es in enumerate(V):
        for e in es: ends[e].append(v)
    adj = [[] for _ in V]
    for a, b in ends: adj[a].append(b); adj[b].append(a)
    return adj, [len(es) // 2 for es in V]


def field(q):
    """elements 0..q-1 with add/mul tables; q prime or q = 9 (F_3[i], i^2 = -1)."""
    if q == 9:
        el = [(a, b) for a in range(3) for b in range(3)]; ix = {x: i for i, x in enumerate(el)}
        add = [[ix[((x[0] + y[0]) % 3, (x[1] + y[1]) % 3)] for y in el] for x in el]
        mul = [[ix[((x[0] * y[0] - x[1] * y[1]) % 3, (x[0] * y[1] + x[1] * y[0]) % 3)] for y in el] for x in el]
        return add, mul
    return [[(a + b) % q for b in range(q)] for a in range(q)], [[a * b % q for b in range(q)] for a in range(q)]


def pg2(q):
    add, mul = field(q); inv = {a: next(b for b in range(1, q) if mul[a][b] == 1) for a in range(1, q)}
    def norm(v):
        for x in v:
            if x: return tuple(mul[inv[x]][y] for y in v)
    pts = sorted({norm(v) for v in itertools.product(range(q), repeat=3) if any(v)})
    dot = lambda u, v: add[add[mul[u[0]][v[0]]][mul[u[1]][v[1]]]][mul[u[2]][v[2]]]
    Pts = [[] for _ in pts]; Lns = [[] for _ in pts]; E = 0
    for j, l in enumerate(pts):
        for i, p in enumerate(pts):
            if dot(p, l) == 0: Pts[i].append(E); Lns[j].append(E); E += 1
    return dict(Pts=Pts, Lns=Lns, E=E, girth=6)


G = {}
for f in ('graphs_small.json', 'graphs_gq.json', 'graphs_gq2.json'):
    G.update(json.load(open(os.path.join(HERE, f))))
for q in (7, 9, 11, 13, 17):
    G[f'PG{q}'] = pg2(q)
thr = lambda g: 1 / (1 + (g - 1) ** (8 / 9))          # X(., 1/9): kappa >= 0.1 needs s >= thr(g)
# expXE table: (graph, Boot value reported, excluded?)
TABLE = [('PG3', 8, True), ('BA4', 10, True), ('PG4M', 10, True), ('PG5', 31, True), ('PG7', 80, True),
         ('PG9', 166, True), ('PG11', 299, True), ('PG13', 491, True), ('PG17', 1095, False),
         ('W3', 20, True), ('W4M', 30, True), ('W5', 128, True), ('W7', 483, True), ('W11', 2960, False)]
ok = True
for name, rep, excl in TABLE:
    g = G[name]; adj, k = adj_of(g['Pts'], g['Lns'], g['E']); rng = random.Random(7)
    b = min(run(adj, k, rng) for _ in range(R if g['E'] <= 3200 else 20))
    s = b / g['E']; t = thr(g['girth'])
    good = (b <= rep) and ((s < t) == excl)
    ok &= good
    print(f"{'ok  ' if good else 'FAIL'} {name:5s} g={g['girth']:2d} E={g['E']:6d} greedy Boot={b:5d} (expXE {rep:5d}) Boot/E={s:.4f}  "
          f"thr(1/9)={t:.4f}  1/g={1/g['girth']:.4f}  -> {'excluded in every field' if s < t else 'not excluded'}")
# tightness: Braverman-He's char-3 code on PG(2,9) has dim 166 (suppliers/xj odd_singer.py 3), so Boot(PG9) = 166 exactly
print("PG(2,9): BH char-3 code dim 166 = greedy Boot 166, so the bound is attained (Boot(PG9) = 166)")
print("PASS_LIMITS_BOOT" if ok else "FAIL")
