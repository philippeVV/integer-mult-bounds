"""A2-triple lifts (S = {s1,s2,s3}, r_s3 = +-r_s1 +- r_s2, F empty): the c_s need not be concurrent (1 extra parameter),
so these are NOT kappa-family members. Exact random members (concurrent and non-concurrent), exact deficient lines on both
sides, and the exact-rho gauge bound mincost (maxdef.py) at these generic points. Usage: a2triple.py seed n"""
import sys, random, json
import numpy as np
from kfam import RT, dot, check, deficient, Fr
from maxdef import mincost
from e8 import *
orb = json.load(open('orbits.json')); S = next(o['rep'] for o in orb['3'] if o['rankS'] == 2)
rnd = random.Random(int(sys.argv[1]))
def rq(): return Fr(rnd.choice([-1, 1]) * rnd.randint(1, 9), rnd.randint(1, 9))
P = [RT[S[0]], RT[S[1]]]
for trial in range(int(sys.argv[2])):
    conc = trial % 2 == 0
    c0 = [rq() for _ in range(8)]
    mu = {s: (rq(), rq(), rq()) for s in S}
    if conc: c = {s: [x + mu[s][0] * y for x, y in zip(c0, RT[s])] for s in S}
    else: c = {s: [x + mu[s][1] * P[0][i] + mu[s][2] * P[1][i] for i, x in enumerate(c0)] for s in S}
    xi = {s: rq() for s in S}; beta = {s: rq() for s in S}
    X = []; Y = []
    for a in range(NL):
        if a in S:
            X.append([xi[a] * v for v in c[a]] + [xi[a]])
            eta = 1 / xi[a] - beta[a] * dot(c[a], RT[a])
            Y.append([beta[a] * v for v in RT[a]] + [eta])
        else:
            X.append(list(RT[a]) + [Fr(0)])
            vals = {-dot(c[s], RT[a]) / 2 for s in S if IP[a][s] == 0}
            assert len(vals) == 1, 'eta inconsistent'
            Y.append([v / 2 for v in RT[a]] + [vals.pop()])
    assert check(X, Y)
    dX, kx = deficient(X); dY, ky = deficient(Y)
    def L(d): return [(np.array([float(x) for x in l]) / np.linalg.norm([float(x) for x in l]), h) for l, h in d.items() if l[0] != 'dim']
    XL, YL = L(dX), L(dY)
    mc = mincost(XL, YL) if kx <= 1 and ky <= 1 else None
    rY = np.linalg.matrix_rank(np.array([[float(x) for x in Y[a]] for a in range(NL) if a not in S]))
    print('concurrent' if conc else 'NON-concurrent', 'rank Y_T', rY, 'X lines', sorted(len(h) for h in dX.values()),
          'Y lines', sorted(len(h) for h in dY.values()), 'kerdim', kx, ky, 'mincost', mc[0] if mc else None, flush=True)
