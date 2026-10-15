"""|S| = 8 (and |S| = 7 with F empty): A = S, B subset of S (F empty), all zero patterns of al/alp on A cap B up to
permutation of S (the relaxed model only sees labels of S, which are orthonormal -> S_k-symmetric; spot-checked below)."""
import sys, random, itertools
from relax import maxcredit
from cliques_lib import random_clique, common
from e8 import *
k = int(sys.argv[1]); rnd = random.Random(k)
while True:
    S = random_clique(k, rnd)
    if len(common(S)) == k: break
print('S', S, 'C(S) = S:', common(S) == sorted(S))
worst = (-1, None); n = 0
for b in range(0, k + 1):
    B = S[:b]
    for n1 in range(b + 1):
        for n2 in range(b + 1 - n1):
            for n3 in range(b + 1 - n1 - n2):
                ZX = B[:n1] + B[n1 + n2:n1 + n2 + n3]; ZY = B[n1:n1 + n2 + n3]
                v, conf = maxcredit(S, B, ZX, ZY, want=True); n += 1
                if v > worst[0]: worst = (v, (b, n1, n2, n3, conf))
                if v > 9: print('CREDIT > 9', b, n1, n2, n3, v, conf, flush=True)
print('cases', n, 'max credit', worst[0], 'at', worst[1])
# symmetry spot check: random full patterns on random B
mx = 0
for t in range(60):
    B = [s for s in S if rnd.random() < 0.6]
    ZX = [s for s in B if rnd.random() < 0.3]; ZY = [s for s in B if rnd.random() < 0.3]
    mx = max(mx, maxcredit(S, B, ZX, ZY))
print('random-pattern spot check max credit', mx)
