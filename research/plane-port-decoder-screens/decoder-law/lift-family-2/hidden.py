"""Are the 'rank-9-on-T' completions 2_1, 5_1, 7_1 lifts with a larger (hidden) S? Find all sparse vectors in col(X):
z with supp(X z) small, via: for each detected-S superset S' = S + extra lines (|S'| <= 8), rank X[V minus S'] == 8."""
import sys, itertools, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from e8 import *
D = np.load('./sols_0_20.npz')
REQ = np.array([[a == b or IP[a][b] == 0 for b in range(NL)] for a in range(NL)])
for k in sys.argv[1:]:
    V, W = D[k]; N = V @ W.T
    dev = (~REQ) & (abs(abs(N * N.T) - 0.25) > 1e-6); deg = dev.sum(1)
    S = [a for a in range(NL) if deg[a] >= 50]
    print(k, 'detected S', S, 'deg of S', [int(deg[a]) for a in S], 'other nonzero deg', sorted(int(d) for a, d in enumerate(deg) if a not in S and d)[-12:])
    for side, M in (('X', V), ('Y', W)):
        found = []
        for r in (0, 1, 2):
            for ext in itertools.combinations([a for a in range(NL) if a not in S], r):
                Sp = S + list(ext)
                T = [a for a in range(NL) if a not in Sp]
                if np.linalg.matrix_rank(M[T], 1e-7) == 8: found.append(ext)
            if found: break
        print('  ', side, 'minimal extensions making rank 8 on complement:', found[:10], len(found))
