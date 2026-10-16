"""rank of the LZD_n Gram W (exact, sympy-free: fraction-free integer elimination via numpy float is fine at this size,
cross-checked by the isotropic/canonical vector k=(3;1^n): all roots lie in k^perp)."""
import itertools, numpy as np
for n in range(8, 16):
    vecs = []
    for T in itertools.combinations(range(1, n + 1), 3):
        v = [0] * (n + 1); v[0] = 1
        for i in T: v[i] = 1
        vecs.append(v)
    for i, j in itertools.combinations(range(1, n + 1), 2):
        v = [0] * (n + 1); v[i] = 1; v[j] = -1; vecs.append(v)
    X = np.array(vecs); eps = np.array([-1] + [1] * n); W = (X * eps) @ X.T
    k = np.array([3] + [1] * n)
    ev = np.linalg.eigvalsh(W.astype(float))
    if n >= 10: assert np.linalg.matrix_rank(W.astype(float)) == n and (ev > 1e-9).sum() == n - 1 and (ev < -1e-9).sum() == 1
    print(f"n={n}: rank W = {np.linalg.matrix_rank(W.astype(float))}, all roots in k^perp: {bool(np.all((X*eps)@k==0))}, "
          f"k.k = {-9+n}, signature (+{int((ev>1e-9).sum())}, -{int((ev<-1e-9).sum())})")
print('PASS_EXPXA_RANKW')
