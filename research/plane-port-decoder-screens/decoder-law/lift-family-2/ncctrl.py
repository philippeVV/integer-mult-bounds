"""Exact control for the non-clique support A = {s, t} (s perp t), B = C2 (24 lines non-orth to both): generic kappa' and the
special level-constant kappa' (pencil line L_s deficient). Prediction: X lines = {e9} (Def dim 7); Y lines = e9' (dim 2)
[+ L_s with Def {486}]; mincost >= 72.  Usage: ncctrl.py seed n"""
import sys, random, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from kfam import RT, dot, check, deficient, Fr
from maxdef import mincost
from e8 import *
rnd = random.Random(int(sys.argv[1]))
def rq(): return Fr(rnd.choice([-1, 1]) * rnd.randint(1, 20), rnd.randint(1, 20))
s = 0; t = next(b for b in range(NL) if IP[s][b] == 0 and b != s)
C2 = [b for b in range(NL) if b not in (s, t) and IP[b][s] != 0 and IP[b][t] != 0]
ok_all = True
for tr in range(int(sys.argv[2])):
    special = tr % 2 == 1
    ka = {s: rq(), t: rq()}; c = rq()
    X = []; Y = []
    for a in range(NL):
        k = ka.get(a, Fr(0))
        kp = (c * Fr(1, 2) * dot(RT[a], RT[s]) if special else rq()) if a in C2 else Fr(0)
        X.append(list(RT[a]) + [k]); Y.append([(1 - k * kp) / 2 * x for x in RT[a]] + [kp])
    assert check(X, Y)
    dX, kx = deficient(X); dY, ky = deficient(Y)
    tov = lambda g: (lambda v: v / np.linalg.norm(v))(np.array([float(x) for x in g]))
    mc = mincost([(tov(g), v) for g, v in dX.items()], [(tov(g), v) for g, v in dY.items()])[0]
    xs = sorted(len(v) for v in dX.values()); ys = sorted(len(v) for v in dY.values())
    ok = kx <= 1 and ky <= 1 and xs == [127] and ys == ([1, 3] if special else [3]) and mc >= 72
    ok_all &= ok
    print(tr, 'special' if special else 'generic', 'X Def sizes', xs, 'Y Def sizes', ys, 'Y classes', sorted(dY.values()), 'mincost', mc, 'OK' if ok else 'FAIL', flush=True)
print('PASS' if ok_all else 'FAIL')
