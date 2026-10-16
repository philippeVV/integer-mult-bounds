"""Exact (mod split prime) sweep of structured Singer-equivariant families on PG(2,17)."""
import sys, json, random, math, time, numpy as np
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from xm_core import *
from geom import singer_coords
_, _, ST = singer_coords(q)
which = int(sys.argv[1]) if len(sys.argv) > 1 else 0
ctx = Ctx(PRIMES[which]); P = ctx.P
z16 = root_of_unity(P, 16, ctx.g); z17 = root_of_unity(P, 17, ctx.g); z3 = root_of_unity(P, 3, ctx.g)
zN = ctx.z
def zp(e): return pow(zN, e % N, P)
eta = lambda x: (zp(x) + zp(17*x) + zp(289*x)) % P      # Gaussian period: generator of the fixed field of σ_17
best = (0, None); rows = []
def rec(name, A, B, mds=False):
    global best
    b = blocks(ctx, A, B); S = int(b.sum())
    r = dict(fam=name, S=S, SE=round(S/E, 4), maxblk=int(b.max()), nz=int((b > 0).sum()), rkA=rank(ctx, A), rkB=rank(ctx, B))
    if mds or S > 150: r['badA'] = mds_bad(ctx, A); r['badB'] = mds_bad(ctx, B)
    rows.append(r); print(json.dumps(r), flush=True)
    if S > best[0]: best = (S, name)
    return S
def grs(xs, vs, kk=k):
    return np.array([[vs[i] * pow(int(xs[i]), r, P) % P for i in range(n)] for r in range(kk)], dtype=np.int64)
def hom_grs(pts, vs, kk=k):   # pts (s,t) -> rows s^{k-1-r} t^r
    return np.array([[vs[i] * pow(int(pts[i][0]), kk-1-r, P) * pow(int(pts[i][1]), r, P) % P for i in range(n)] for r in range(kk)], dtype=np.int64)
rng = random.Random(1)
ones = [1]*n
# --- F1: frequency codes (Chebotarev -> S = 81) -------------------------------------------------
for t in range(3):
    FA = rng.sample(range(N), k); FB = rng.sample(range(N), k)
    A = np.array([[zp(f*d) for d in D] for f in FA]); B = np.array([[zp(f*d) for d in D] for f in FB]); rec(f'freq_rand{t}', A, B)
# --- F2: Galois-twisted GRS at roots of unity x_d = ζ_N^{a d}, weights in the σ_17-fixed field --------
units = [u for u in range(1, N)]
orbreps = sorted({min(u*17**j % N for j in range(3)) for u in units})   # Z_N^*/<17>, 102 reps
for wname, wfun in [('w1', lambda d: 1), ('eta', lambda d: eta(d)), ('eta2', lambda d: eta(d*d)), ('eta_d+1', lambda d: eta(d) + 1)]:
    v = [wfun(d) % P for d in D]
    if 0 in v: continue
    A = grs([zp(d) for d in D], v)
    for u in [1, 17, 289, -1, 2, 3, 5] + rng.sample(orbreps, 6):
        B = grs([zp(u*d) for d in D], [eta(u*d) if wname.startswith('eta') else 1 for d in D]) if wname != 'w1' else grs([zp(u*d) for d in D], ones)
        if wname == 'eta2': B = grs([zp(u*d) for d in D], [eta(u*u*d*d) for d in D])
        if wname == 'eta_d+1': B = grs([zp(u*d) for d in D], [(eta(u*d)+1) % P for d in D])
        rec(f'GRS_rootsN_{wname}_sigma{u}', A, B)
# --- F3: GRS on the pencil P^1(F_17), lifted to char 0 -------------------------------------------
pts = [ST[d] for d in D]
lg = {}; g17 = 3
for e in range(16): lg[pow(g17, e, 17)] = e
def teich(x, u=1): return 0 if x == 0 else pow(z16, (lg[x]*u) % 16, P)
def addlift(x, u=1): return pow(z17, (x*u) % 17, P)
for u in [1, 3, 5, 7, 9, 11, 13, 15]:
    A = hom_grs([(teich(s), teich(t)) for s, t in pts], ones)
    B = hom_grs([(teich(s, u), teich(t, u)) for s, t in pts], ones)
    rec(f'teich_homGRS_twist{u}', A, B)
for u in range(1, 17):
    # additive-character lift of the affine coordinate; ∞ -> 0
    xs = [0 if t == 0 else addlift(s) for s, t in pts]
    xsu = [0 if t == 0 else addlift(s, u) for s, t in pts]
    rec(f'addlift_GRS_sigma{u}', grs(xs, ones), grs(xsu, ones))
# --- F4: Gaussian-period (σ_17-fixed-field) Cauchy/Vandermonde-type codes ------------------------
for t in range(4):
    a = rng.sample(orbreps, k); c = rng.sample(range(N), n)
    A = np.array([[eta(ai * (d + c[j]) ) for j, d in enumerate(D)] for ai in a])
    for u in [1, -1, 2]:
        B = np.array([[eta(u*ai * (d + c[j])) for j, d in enumerate(D)] for ai in a])
        rec(f'eta_mat{t}_sigma{u}', A, B)
# --- F5: semilinear Frobenius analogue: B = π∘σ_u(A) with π: d -> 17 d (multiplier), A over Q(ζ_N) ----
perm = [D.index((17*d) % N) for d in D]
for t in range(3):
    A = np.array([[ (zp(f*d) + 2*zp(f*d*d)) % P for d in D] for f in rng.sample(range(1, N), k)])
    for u in [1, -1, 2, 17]:
        sig = np.array([[0]*n]*k)
        # apply σ_u entrywise is only meaningful if we know the exponent form; rebuild
    FA = rng.sample(range(1, N), k)
    for u in [-1, 2, 3]:
        A = np.array([[(zp(f*d) + zp(f*d*d)) % P for d in D] for f in FA])
        B = np.array([[(zp(u*f*d) + zp(u*f*d*d)) % P for d in D] for f in FA])[:, perm]
        rec(f'semilin_quad{t}_sigma{u}_pi', A, B)
        rec(f'semilin_quad{t}_sigma{u}', A, np.array([[(zp(u*f*d) + zp(u*f*d*d)) % P for d in D] for f in FA]))
# --- F6: generic baselines ------------------------------------------------------------------------
for t in range(2):
    A = np.array([[rng.randrange(P) for _ in range(n)] for _ in range(k)]); B = np.array([[rng.randrange(P) for _ in range(n)] for _ in range(k)])
    rec(f'generic{t}', A, B)
# conjugate pair
A = np.array([[rng.randrange(P) for _ in range(n)] for _ in range(k)]); rec('A=B generic', A, A)
print('BEST', best, 'threshold M6*E=1005, Boot=1095')
json.dump(rows, open(f'fam_P{which}.json', 'w'))
assert best[0] == 81 and best[0] < 1005, best   # limits-1 replay: every exact equivariant family <= k^2 = 81 << M_6 E = 1005
print('PASS_EXPXM_FAM')
