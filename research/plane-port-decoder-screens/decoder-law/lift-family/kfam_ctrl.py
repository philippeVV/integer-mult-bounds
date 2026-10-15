"""Exact control: random kappa-family members are exact masked completions; their deficient lines are exactly the predicted
ones: X: e9 (Def = span(q_A)^perp) and g_f = (ka_f r_f, -2 al_f) for f in A with (A cap orth f = {} or al_f = 0);
Y symmetric with B. Usage: kfam_ctrl.py seed"""
import sys, random
from kfam import *
from cliques_lib import random_clique
rnd = random.Random(int(sys.argv[1]))
def span_perp(labs):
    return [h for h in range(1, 512) if all(dotq(l, h) == 0 for l in labs)]
def predict(A, P, side):
    L = {}
    e9 = tuple([Fr(0)] * 8 + [Fr(1)])
    D = [h for h in span_perp([LAB[a] for a in A]) if h != 511]
    if D: L[e9] = sorted(D)
    for f in A:
        x, k, y, kp = P[f]
        al, ka = (x, k) if side == 'X' else (y, kp)
        if all(IP[f][a] != 0 for a in A if a != f) or al == 0:
            g = normline([ka * c for c in RT[f]] + [-2 * al])
            L.setdefault(g, []).append(LAB[f]); L[g] = sorted(L[g])
    return L
tests = 0
for trial in range(int(sys.argv[2])):
    k = rnd.choice([1, 2, 3, 4, 5, 6, 7, 8, 8, 8])
    S = random_clique(k, rnd)
    from cliques_lib import common
    C = common(S)
    mode = rnd.choice(['subS', 'subC', 'small'])
    if mode == 'subS': B = [s for s in S if rnd.random() < 0.7]
    elif mode == 'subC': B = [c for c in C if rnd.random() < 0.5][:8]
    else: B = rnd.sample(C, min(len(C), rnd.randint(0, 3)))
    AB = [a for a in S if a in B]
    al = {a: Fr(0) for a in AB if rnd.random() < 0.3}
    alp = {a: Fr(0) for a in AB if a not in al and rnd.random() < 0.3}
    X, Y, P = member(S, B, rnd, al, alp)
    assert check(X, Y)
    from orbits import clean
    cl = not clean(sum(1 << a for a in S)) and not clean(sum(1 << b for b in B))
    # pairing identities: e9.e9' = 1, e9.h'_f = -2 alp_f, g_f.e9' = -2 al_f, g_f.h'_f = 2 (ka_f kap_f + 2 al_f alp_f) = 2
    for f in set(S) & set(B):
        x, k_, y, kp = P[f]
        g = [k_ * c for c in RT[f]] + [-2 * x]; hp = [kp * c for c in RT[f]] + [-2 * y]
        assert dot(g, hp) == 2 and dot(g, [0]*8 + [1]) == -2 * x and dot([0]*8 + [1], hp) == -2 * y
    dX, kx = deficient(X); dY, ky = deficient(Y)
    pX = predict(S, P, 'X'); pY = predict(B, P, 'Y')
    okX = {l: sorted(h) for l, h in dX.items()} == pX
    okY = {l: sorted(h) for l, h in dY.items()} == pY
    print(trial, '|S|', k, 'B', len(B), '|C|', len(C), 'zeros', len(al), len(alp), 'kerdim', kx, ky, 'clean', cl, 'X ok', okX, 'Y ok', okY, flush=True)
    if cl: assert okX and okY, (dX, pX, dY, pY)
    else: print('   (B not Clean: prediction not claimed)')
    tests += 1
print('PASS', tests)
