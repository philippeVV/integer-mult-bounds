"""Exact members of the non-concurrent A2-triple family (normal form of a2lib.py) incl. special parameters
(beta_s = 0; merged h'_i = h'_j = L_ij; mu choices making e9.h'_i or g_i.h'_j vanish): exact deficient lines on both sides
vs the prediction of a2credit.py, and the exact-rho gauge floor mincost >= 72.  Usage: a2ctrl.py seed n"""
import sys, random, numpy as np
from a2lib import *
from kfam import check, deficient
from maxdef import mincost
rnd = random.Random(int(sys.argv[1]))
def rq(): return Fr(rnd.choice([-1, 1]) * rnd.randint(1, 12), rnd.randint(1, 12))
K = ylines_T()
L = {}
for h, k in K.items():
    if len(k) == 1: L[tuple(sorted(set(CLASS[h]) & set(S)))] = k[0]
ok_all = True
for tr in range(int(sys.argv[2])):
    kind = tr % 4
    mu = {x: rq() for x in S}; beta = {x: rq() for x in S}; eta = {x: rq() for x in S}
    merged = None
    if kind == 1:
        x = rnd.choice(S); beta[x] = Fr(0)
    if kind in (2, 3):
        pr = rnd.choice(list(L)); merged = pr; Lv = L[pr]
        for x in pr: eta[x] = -beta[x] * dot(RT[x], Lv[:8]) / Lv[8]
        if kind == 3: mu = {x: Fr(0) for x in S}
    mkc = lambda: {s1: [mu[s1] * v for v in RT[s1]], s2: [mu[s2] * v for v in RT[s2]], s3: [mu[s3] * v + w for v, w in zip(RT[s3], n)]}
    c = mkc()
    while any(beta[a] * dot(c[a], RT[a]) + eta[a] == 0 for a in S):     # N_ss = 1 needs den != 0: re-draw mu (keeps the merge)
        mu = {x: rq() for x in S}; c = mkc()
    X = []; Y = []
    for a in range(NL):
        if a in S:
            den = beta[a] * dot(c[a], RT[a]) + eta[a]
            xi = 1 / den
            X.append([xi * v for v in c[a]] + [xi]); Y.append([beta[a] * v for v in RT[a]] + [eta[a]])
        else:
            X.append(list(RT[a]) + [Fr(0)]); Y.append(YT[a])
    assert check(X, Y)
    dX, kx = deficient(X); dY, ky = deficient(Y)
    xs = sorted(len(h) for h in dX.values()); ys = sorted(len(h) for h in dY.values())
    pred_y = [1, 3] if merged and all(True for _ in [0]) else [1, 1, 1]
    tov = lambda g: (lambda v: v / np.linalg.norm(v))(np.array([float(x) for x in g]))
    mc = mincost([(tov(g), v) for g, v in dX.items()], [(tov(g), v) for g, v in dY.items()])[0]
    ok = kx <= 1 and ky <= 1 and xs == [1, 1, 1, 63] and ys == pred_y and mc >= 72
    ok_all &= ok
    print(tr, 'kind', kind, 'merged', merged, 'X Def sizes', xs, 'Y Def sizes', ys, 'pred', pred_y, 'mincost', mc, 'OK' if ok else 'FAIL', flush=True)
print('PASS' if ok_all else 'FAIL')
