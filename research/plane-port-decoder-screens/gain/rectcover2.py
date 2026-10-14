"""Frame (G/2) + rectangle corrections, EXACT. Planes: frame-compatible (keeps frame terms at cost 8) found by local
search, or loaded. Lost entries (odd pair, planes not orthogonal) need value -G_ab/2. Correction term = rectangle
R x A of uncovered lost entries with rank-1 sign pattern, cost = F2 pairing rank of the plane spans. Greedy partition
(ratio entries/cost), then exact verification with Fractions. usage: rectcover.py seed mode(frame|free|file:path) maxcost"""
import sys, json, random, os
from fractions import Fraction as Fr
from base import *
HERE = os.path.dirname(os.path.abspath(__file__))
emap = {int(k): v for k, v in json.load(open(os.path.join(HERE, '../lower/e8map.json'))).items()}
ip = lambda a, b: sum(x * y for x, y in zip(a, b)) // 4
Rt = [emap[q] for q in P]; frame = [emap[x] for x in sorted(set(emap) - set(P))]
c = [[ip(Rt[a], f) for f in frame] for a in range(n)]
D = [[a for a in range(n) if c[a][k]] for k in range(8)]
G = [[ip(Rt[a], Rt[b]) for b in range(n)] for a in range(n)]
seed = int(sys.argv[1]); mode = sys.argv[2]; maxcost = int(sys.argv[3]); rng = random.Random(seed)
def span(vs):
    out = [0]
    for b in reduce(vs): out += [v ^ b for v in out]
    return set(out)
if mode.startswith('file:'):
    d = json.load(open(mode[5:])); d = d[0] if isinstance(d, list) else d; r = d.get('r') or d.get('R')
else:
    W = [span([P[a] for a in D[k]]) for k in range(8)]
    dom = [[x for x in rchoices(P[a]) if mode == 'free' or all(x in W[k] for k in range(8) if c[a][k])] for a in range(n)]
    r = [rng.choice(dom[a]) for a in range(n)]
    def orth(i, ri, j, rj): return dot(P[i], rj) == 0 and dot(ri, P[j]) == 0 and dot(ri, rj) == 0
    def fdeg(i, ri): return sum(orth(i, ri, j, r[j]) for j in NB[i])
    import math; T0 = 2.0
    for it in range(60000):
        i = rng.randrange(n); ri = rng.choice(dom[i]); dl = fdeg(i, ri) - fdeg(i, r[i])
        if dl >= 0 or rng.random() < math.exp(dl / T0): r[i] = ri
        T0 = max(0.05, T0 * 0.9999)
F = free_matrix(r); E = planes(r)
lost = {(a, b) for a in range(n) for b in range(n) if a != b and G[a][b] % 2 and not F[a][b]}
fcost = [term_cost(r, D[k], D[k]) for k in range(8)]
print('seed', seed, mode, 'free pairs', sum(map(sum, F)) // 2, 'lost entries', len(lost), 'frame costs', fcost, flush=True)
unc = set(lost)
rows_of = {}
def cost(R, A): return pair_rank([w for i in R for w in E[i]], [w for j in A for w in E[j]])
def grow(a, b):
    R, A = [a], [b]; sig = {a: 1}; tau = {b: G[a][b]}; cst = cost(R, A)
    if cst > maxcost: return None
    while True:
        # candidate rows / cols
        cr = [x for x in range(n) if x not in sig and all((x, y) in unc or F[x][y] for y in A)
              and len({G[x][y] * tau[y] for y in A if (x, y) in unc}) == 1]
        cc = [y for y in range(n) if y not in tau and all((x, y) in unc or F[x][y] for x in R)
              and len({G[x][y] * sig[x] for x in R if (x, y) in unc}) == 1]
        rng.shuffle(cr); rng.shuffle(cc)
        added = False
        opts = [('r', x) for x in cr] + [('c', y) for y in cc]
        rng.shuffle(opts)
        for kind, z in opts:
            R2, A2 = (R + [z], A) if kind == 'r' else (R, A + [z])
            c2 = cost(R2, A2)
            if c2 <= cst:
                if kind == 'r':
                    y0 = next(y for y in A if (z, y) in unc); sig[z] = G[z][y0] * tau[y0]; R = R2
                else:
                    x0 = next(x for x in R if (x, z) in unc); tau[z] = G[x0][z] * sig[x0]; A = A2
                added = True; break
        if not added: return R, A, cst, sig, tau
terms = []
while unc:
    best = None
    pool = rng.sample(sorted(unc), min(40, len(unc)))
    for a, b in pool:
        g = grow(a, b)
        if g is None: continue
        R, A, cst, sig, tau = g
        sc = sum((x, y) in unc for x in R for y in A) / max(cst, 1)
        if best is None or sc > best[0]: best = (sc, g)
    if best is None:  # fall back: single entry
        a, b = next(iter(unc)); g = ([a], [b], cost([a], [b]), {a: 1}, {b: G[a][b]}); best = (0, g)
    R, A, cst, sig, tau = best[1]
    for x in R:
        for y in A: unc.discard((x, y))
    terms.append(({x: Fr(-sig[x], 2) for x in R}, {y: Fr(tau[y]) for y in A}))
dec = [({a: Fr(c[a][k], 4) for a in D[k]}, {a: Fr(c[a][k]) for a in D[k]}) for k in range(8)] + terms
bad, C, costs = verify(r, dec)
from collections import Counter
print('corr terms', len(terms), 'corr cost hist', sorted(Counter(costs[8:]).items()), 'sizes', sorted(Counter(len(d)*len(t) for d,t in terms).items())[-5:])
print('EXACT bad', bad, 'C', C, flush=True)
json.dump(dict(r=r, C=C, bad=bad, costs=costs, terms=[[{str(k): str(v) for k, v in d.items()}, {str(k): str(v) for k, v in t.items()}] for d, t in dec]),
          open(f'rectF_{seed}_{mode.split(":")[0]}_{maxcost}.json', 'w'))
