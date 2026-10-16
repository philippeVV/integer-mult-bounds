import sys, numpy as np, itertools
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from xg_core import *
V, nP, E, g = load('PG3')
for f in sys.argv[1:]:
    Qs = [orth(np.asarray(Q, complex)) for Q in np.load(f, allow_pickle=True)]
    for side in (0, 1):
        ld, err = flat_gauge(V, nP, E, Qs, side, 5000, 1e-14)
        n1, w1 = nullity(V, nP, E, Qs, logd=ld)
        Ys = [orth(np.exp(ld[es])[:, None] * Q) for es, Q in zip(V, Qs)]
        mg = min(abs(np.linalg.det(Y[list(S)])) for Y in Ys for S in itertools.combinations(range(4), 2))
        print(f, 'flat side', side, 'eig', ['%.2e' % x for x in w1[:11]], 'min 2x2 minor %.2e' % mg)
