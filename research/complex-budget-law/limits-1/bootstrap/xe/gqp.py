"""W(q), q prime: points of PG(3,q), totally isotropic lines; and bootstrap-ready JSON (graphs_gq.json)."""
import sys, json, itertools, os
out = {}
for q in map(int, sys.argv[1:]):
    def norm(v):
        for x in v:
            if x % q:
                inv = pow(x, q - 2, q); return tuple(y * inv % q for y in v)
    pts = sorted({norm(v) for v in itertools.product(range(q), repeat=4) if any(v)})
    idx = {p: i for i, p in enumerate(pts)}
    B = lambda x, y: (x[0]*y[1] - x[1]*y[0] + x[2]*y[3] - x[3]*y[2]) % q
    lines = {}
    for i, P in enumerate(pts):
        for Q in pts[i + 1:]:
            if B(P, Q): continue
            key = tuple(sorted({idx[norm(tuple((a*P[t] + b*Q[t]) % q for t in range(4)))] for a in range(q) for b in range(q) if a or b}))
            lines.setdefault(key, len(lines))
    Ls = sorted(lines, key=lines.get); Pinc = [[] for _ in pts]
    for l, s in enumerate(Ls):
        for p in s: Pinc[p].append(l)
    eid = {}; Pts = []
    for v, ls in enumerate(Pinc):
        Pts.append([eid.setdefault((v, l), len(eid)) for l in ls])
    Lns = [[eid[(p, l)] for p in s] for l, s in enumerate(Ls)]
    out[f'W{q}'] = dict(Pts=Pts, Lns=Lns, E=len(eid), girth=8); print(q, len(pts), len(Ls), len(eid), flush=True)
json.dump(out, open(os.environ.get('OUT', 'graphs_gq.json'), 'w'))
