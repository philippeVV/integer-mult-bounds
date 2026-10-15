"""P120 searches for decoders with C <= 71 built from cheap label-subspace supports.
Hyperplane classes on P120 (h in F2^9): (even h, 64 lines, rho 7) x120 [E7+A1 = perp(f) + f],
(even, 56, 7) x135 [D8 type], (odd, 56, 8) x120 [non-orth to f], (odd, 64, 8) x135."""
import sys, json, time, itertools
import numpy as np
from lmg import *
from pysat.solvers import Cadical153
from pysat.formula import IDPool
idx, n, REQ = setup(FULL)
HYP = {h: sum(1 << a for a in range(NL) if dotq(LAB[a], h) == 0) for h in range(1, 512)}
CL = {}
for h, m in HYP.items():
    if m: CL.setdefault((bin(h).count('1') % 2, bin(m).count('1'), rho(m, m)), []).append(h)
def masks(sup):
    K = len(sup); Vm = np.zeros((n, K), bool); Wm = np.zeros((n, K), bool)
    for k, (R, A) in enumerate(sup):
        for a in bits(R): Vm[a, k] = True
        for a in bits(A): Wm[a, k] = True
    return Vm, Wm
def boolsat(sup, sym=False):
    K = len(sup); pool = IDPool(); cl = []
    X = {(a, k): pool.id(('x', a, k)) for k in range(K) for a in bits(sup[k][0])}
    Y = X if sym else {(b, k): pool.id(('y', b, k)) for k in range(K) for b in bits(sup[k][1])}
    for a in range(n):
        for b in range(n):
            if not REQ[a, b]: continue
            ks = [k for k in range(K) if (a, k) in X and (b, k) in Y]
            if a == b:
                if not ks: return False
                ps = []
                for k in ks:
                    p = pool.id(('p', a, b, k)); ps.append(p); cl += [[-p, X[a, k]], [-p, Y[b, k]]]
                cl.append(ps)
            elif ks:
                if len(ks) == 1: cl.append([-X[a, ks[0]], -Y[b, ks[0]]]); continue
                ps = []
                for k in ks:
                    p = pool.id(('p', a, b, k)); ps.append(p)
                    cl += [[-p, X[a, k]], [-p, Y[b, k]], [p, -X[a, k], -Y[b, k]]]
                for i, p in enumerate(ps): cl.append([-p] + [q for j, q in enumerate(ps) if j != i])
    s = Cadical153(bootstrap_with=cl); r = s.solve(); s.delete(); return r
def cost(sup): return sum(rho(R, A) for R, A in sup)
if __name__ == '__main__':
    mode, N, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    rng = np.random.default_rng(seed); t0 = time.time(); stats = dict(tried=0, sat=0, best=None)
    E64 = CL[(0, 64, 7)]; E56 = CL[(0, 56, 7)]; O56 = CL[(1, 56, 8)]; O64 = CL[(1, 64, 8)]
    for it in range(N):
        if mode == 'e64x9': hs = list(rng.choice(E64, 9, replace=False)); sup = [(HYP[h], HYP[h]) for h in hs]
        elif mode == 'e64x10': hs = list(rng.choice(E64, 10, replace=False)); sup = [(HYP[h], HYP[h]) for h in hs]
        elif mode == 'e7mix10':   # 10 rho-7 hyperplanes from both even classes
            hs = list(rng.choice(E64 + E56, 10, replace=False)); sup = [(HYP[h], HYP[h]) for h in hs]
        elif mode == 'asym9':     # 9 asym pairs (U, U') hyperplanes, any classes, cost <= 71
            hs = [(int(rng.choice(E64 + E56 + O56 + O64)), int(rng.choice(E64 + E56 + O56 + O64))) for _ in range(9)]
            sup = [(HYP[a], HYP[b]) for a, b in hs]
        elif mode == 'mix9':      # 9 symmetric hyperplanes, any classes, cost <= 71
            hs = list(rng.choice(E64 + E56 + O56 + O64, 9, replace=False)); sup = [(HYP[h], HYP[h]) for h in hs]
        c = cost(sup)
        if c > 71: continue
        stats['tried'] += 1
        sym = mode != 'asym9'
        if not boolsat(sup, sym=sym): continue
        stats['sat'] += 1
        r = lm(FULL, *masks(sup), seed=it, tries=8, maxfev=300, sym=sym)
        res = r[0] if r else None
        if res is not None and (stats['best'] is None or res < stats['best'][0]): stats['best'] = (res, [int(x) if not isinstance(x, tuple) else x for x in hs], c)
        print(it, 'cost', c, 'res', None if res is None else round(res, 4), flush=True)
        if res is not None and res < 1e-8: print('FOUND', hs, c, flush=True)
    print('stats', stats['tried'], stats['sat'], 'best', stats['best'], 'time', round(time.time() - t0))
