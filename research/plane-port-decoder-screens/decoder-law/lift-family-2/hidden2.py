"""Sign-twisted completions 0_1, 13_1: find sparse vectors z in col(X) (X z supported on few lines) through the deficient
lines, then test the kappa-family normal form with S' = supp(X z)."""
import sys, numpy as np, itertools
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from e8 import *
D = np.load('./sols_0_20.npz')
CLS = {h: np.array([dotq(LAB[a], h) for a in range(NL)], bool) for h in range(1, 511)}
REQ = np.array([[a == b or IP[a][b] == 0 for b in range(NL)] for a in range(NL)])
R = np.array(ROOT, float) / 2
for k in sys.argv[1:]:
    V, W = D[k]
    for side, M in (('X', V), ('Y', W)):
        best = None
        for h, O in CLS.items():
            u, s, vt = np.linalg.svd(M[O], full_matrices=False)
            if s[-1] < 1e-7 * s[0]:
                z = vt[-1]; v = M @ z; sup = np.where(np.abs(v) > 1e-7 * np.abs(v).max())[0]
                if best is None or len(sup) < len(best[1]): best = (z, list(sup))
        S = best[1]
        T = [a for a in range(NL) if a not in S]
        pw = all(IP[a][b] != 0 for a, b in itertools.combinations(S, 2))
        print(k, side, 'sparsest deficient-line support', S, 'size', len(S), 'pairwise odd', pw, 'rank on complement', np.linalg.matrix_rank(M[T], 1e-7))
    # kappa-family fit: torus on T, then X_T = R_T A
    N = V @ W.T
    Np = N.copy()
    # torus from rigid T (expSM/param.py)
    sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__))); from param import gauge_fix
    X2, Y2, d, err, fit = gauge_fix(V, W, S)
    alpha = np.array([X2[a, :8] @ R[a] / 2 for a in range(NL)]); alp = np.array([Y2[a, :8] @ R[a] / 2 for a in range(NL)])
    resX = max(np.abs(X2[a, :8] - alpha[a] * R[a]).max() for a in range(NL)); resY = max(np.abs(Y2[a, :8] - alp[a] * R[a]).max() for a in range(NL))
    print('   torus err on T', f'{err:.1e}', 'X first-8 = alpha r residual', f'{resX:.1e}', 'Y first-8 = alpha\' r residual', f'{resY:.1e}')
    ka = X2[:, 8]; kap = Y2[:, 8]
    print('   supp kappa', list(np.where(np.abs(ka) > 1e-7)[0]), 'supp kappa\'', list(np.where(np.abs(kap) > 1e-7)[0]))
    print('   2 al al\' on S', np.round(2 * alpha[S] * alp[S], 4), ' ka ka\' on S', np.round(ka[S] * kap[S], 4))
