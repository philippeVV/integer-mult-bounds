"""F8-line ports: masked decoder with every term supported inside one F8-plane of F8^3 (projective line of PG(2,8)).
Term cost <= F2 pairing rank of that plane: 6 (nondegenerate) or 3 (degenerate, radical = isotropic F8-line).
Levenberg-Marquardt feasibility (floating point, uncertified) for a chosen multiset of planes.
usage: e_planes.py trial n_deg_copies n_nondeg [restarts] [seed]"""
import sys, itertools, numpy as np
from scipy.optimize import least_squares
from scipy.sparse import coo_matrix
from common import *
import f8
trial, kdeg, knd = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
R = int(sys.argv[4]) if len(sys.argv) > 4 else 2; seed = int(sys.argv[5]) if len(sys.argv) > 5 else 0
F, free = f8.build(1, trial); n = len(F); Fm = np.array(free)
mul, h = f8.mul, f8.h
def proj(u):  # normalise F8 vector to a projective point
    for c in u:
        if c:
            inv = next(x for x in range(1, 8) if mul(c, x) == 1); return tuple(mul(inv, y) for y in u)
pts = sorted({proj(v) for v in itertools.product(range(8), repeat=3) if any(v)})
assert len(pts) == 73
port_pt = [proj(u) for u in F]
planes = []  # projective lines = h-orthogonal complements of points p (h nondegenerate on F8^3)
for p in pts:
    members = [a for a in range(n) if h(F[a], list(p)) == 0]
    iso = h(list(p), list(p)) == 0
    # plane p^perp is degenerate iff p is isotropic (p in p^perp)
    planes.append(dict(normal=p, ports=members, degenerate=iso))
deg = [pl for pl in planes if pl['degenerate']]; nd = [pl for pl in planes if not pl['degenerate']]
# F2 pairing rank of each plane (sanity: 3 for degenerate, 6 otherwise)
def f2vecs(p):
    out = []
    for v in itertools.product(range(8), repeat=3):
        if h(list(v), list(p)) == 0: out.append(v)
    return out
def to_int(v):  # coordinates of F8^3 in the self-dual basis -> 9-bit F2 vector (dot = Tr h)
    x = 0
    for i, c in enumerate(v):
        for j, b in enumerate(f8.SD):
            # coefficient of b in c: Tr(c*b) since basis is self-dual
            if f8.tr(mul(c, b)): x |= 1 << (3 * i + j)
    return x
for pl in deg[:2] + nd[:2]:
    pl['rho'] = pair_rank([to_int(v) for v in f2vecs(pl['normal'])], [to_int(v) for v in f2vecs(pl['normal'])])
print('planes: degenerate', len(deg), 'nondegenerate', len(nd), 'sample rho', [(pl['degenerate'], pl['rho']) for pl in deg[:2] + nd[:2]])
print('ports per degenerate plane', sorted(len(pl['ports']) for pl in deg), ' per nondeg plane (min/max)', min(len(pl['ports']) for pl in nd), max(len(pl['ports']) for pl in nd))
rng0 = np.random.default_rng(seed)
sel = (sorted(range(len(nd)), key=lambda i: -len(nd[i]['ports']))[:knd] if seed == 99 else
       list(rng0.choice(len(nd), size=min(knd, len(nd)), replace=False)))
chosen = [pl for pl in deg for _ in range(kdeg)] + [nd[i] for i in sel]
K = len(chosen); cost = sum(3 if pl['degenerate'] else 6 for pl in chosen)
cover = set(a for pl in chosen for a in pl['ports'])
print('terms', K, 'cost bound', cost, 'ports covered', len(cover), '/', n)
# variables: X[a,k], Y[a,k] for a in plane k
idx = [(a, k) for k, pl in enumerate(chosen) for a in pl['ports']]
pos = {p: i for i, p in enumerate(idx)}; nv = len(idx)
rows = [(a, a) for a in range(n)] + [(a, b) for a in range(n) for b in range(n) if a != b and not Fm[a, b]]
# only rows where some term can contribute matter for the residual of required zeros; all rows kept
ra = np.array([r[0] for r in rows]); rb = np.array([r[1] for r in rows]); tgt = (ra == rb).astype(float)
memb = [set(pl['ports']) for pl in chosen]
pairs = [[] for _ in rows]
for i, (a, b) in enumerate(rows):
    for k in range(K):
        if a in memb[k] and b in memb[k]: pairs[i].append((pos[(a, k)], pos[(b, k)]))
I = np.array([i for i, ps in enumerate(pairs) for _ in ps]); PX = np.array([p for ps in pairs for p, _ in ps]); PY = np.array([q for ps in pairs for _, q in ps])
def fun(z):
    x, y = z[:nv], z[nv:]
    return np.bincount(I, weights=x[PX] * y[PY], minlength=len(rows)) - tgt
def jac(z):
    x, y = z[:nv], z[nv:]
    return coo_matrix((np.concatenate([y[PY], x[PX]]), (np.concatenate([I, I]), np.concatenate([PX, nv + PY]))), shape=(len(rows), 2 * nv)).tocsr()
best = 1e9
for rs in range(R):
    rng = np.random.default_rng(1000 * seed + rs)
    z0 = rng.standard_normal(2 * nv)
    s = least_squares(fun, z0, jac=jac, method='trf', max_nfev=400, xtol=1e-13, ftol=1e-13)
    best = min(best, 2 * s.cost)
    print('restart', rs, 'residual %.3g' % (2 * s.cost), flush=True)
    if best < 1e-12: break
print('RESULT trial', trial, 'kdeg', kdeg, 'knd', knd, 'K', K, 'cost', cost, 'best residual %.3g' % best)
