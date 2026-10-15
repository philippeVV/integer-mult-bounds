import sys, numpy as np, itertools
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from e8 import *
from param import gauge_fix
D = np.load('./sols_0_20.npz'); R = np.array(ROOT, float) / 2
SS = {'0_1': [5, 61, 89, 93, 95, 98, 104], '13_1': [3, 16, 30, 40, 41, 74, 106, 116]}
for k, S in SS.items():
    V, W = D[k]; X2, Y2, d, err, fit = gauge_fix(V, W, S)
    ka = X2[S, 8]; kap = Y2[S, 8]; c = np.sqrt(np.abs(kap).max() / np.abs(ka).max()); ka *= c; kap /= c
    al = np.array([X2[a, :8] @ R[a] / 2 for a in S]); alp = np.array([Y2[a, :8] @ R[a] / 2 for a in S])
    print(k, 'kappa', np.array2string(ka, precision=3), "kappa'", np.array2string(kap, precision=3))
    print('   ka*kap diag', np.array2string(ka * kap, precision=3), ' 2 al al\'', np.array2string(2 * al * alp, precision=4))
    N = V @ W.T
    prods = [N[a, b] * N[b, a] for a, b in itertools.combinations(S, 2)]
    print('   N_ab N_ba on S pairs', np.array2string(np.array(prods), precision=4))
    print('   rank N', np.linalg.matrix_rank(N, 1e-8), ' sv tail', np.linalg.svd(N, compute_uv=False)[7:10])
