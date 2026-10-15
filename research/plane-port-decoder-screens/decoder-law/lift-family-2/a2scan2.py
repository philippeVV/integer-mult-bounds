import collections
from a2lib import *
import numpy as np
K = ylines_T()
E7 = {x: next(h for h in CLASS if set(CLASS[h]) & set(S) == {x} and len(K[h]) == 2) for x in S}
PR = {}
for h, k in K.items():
    if len(k) == 1: PR[tuple(sorted(set(CLASS[h]) & set(S)))] = (h, k[0])
def rk(vs): return np.linalg.matrix_rank(np.array([[float(x) for x in v] for v in vs]), 1e-10)
for x in S:
    print('K_%d' % x, 'class', E7[x], 'label q', LAB[x], E7[x] == LAB[x])
for pr, (h, L) in PR.items():
    print('L', pr, 'class', h, 'in K_i:', {x: rk(K[E7[x]] + [L]) == 2 for x in S}, 'h.q_s', [dotq(h, LAB[x]) for x in S])
print('K_i cap K_j dims', {(x, y): 4 - rk(K[E7[x]] + K[E7[y]]) for x, y in itertools.combinations(S, 2)})
print('span of the 3 L lines rank', rk([L for _, L in PR.values()]), '; L.e9-component (gamma):', [str(L[8]) for _, L in PR.values()])
print('labels q_S independent:', len({0, LAB[s1], LAB[s2], LAB[s3], LAB[s1]^LAB[s2], LAB[s1]^LAB[s3], LAB[s2]^LAB[s3], LAB[s1]^LAB[s2]^LAB[s3]}) == 8)
print('pair-class labels', {pr: h for pr, (h, L) in PR.items()}, 'q', [LAB[x] for x in S], 'q.q pairwise', [dotq(LAB[a], LAB[b]) for a, b in itertools.combinations(S, 2)])
