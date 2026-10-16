"""Second-order model: per-vertex split of the radial derivative at the gauge optimum, r_v = (1/16) z_v^T (H0 z)_v. Is each r_v <= 0?"""
from model import *
from xr_tools import setup
for base in ['sic', '4']:
    S, Q0, W, G, C, D, lmin, rng = setup('W3', base, 0); E = S.E; V = S.V; nP = S.nP
    X, Y = XYmats(V, nP, E, W); I = np.eye(E); H0 = np.block([[4*Y, 3*(I+Y@X)], [3*(I+X@Y), 4*X]])
    L = -np.linalg.pinv(G, rcond=1e-10) @ C; pos = 0; tot = 0
    for k in range(20):
        d = rng.standard_normal(E); x = L @ d; u, v = x, x + d
        z = np.r_[0.5*(I-X)@u, 0.5*(I-Y)@v]; Hz = H0 @ z
        r = np.array([z[es] @ Hz[es] if i < nP else z[E + np.array(es)] @ Hz[E + np.array(es)] for i, es in enumerate(V)]) / 16
        pos += (r > 1e-12).sum(); tot += len(r)
        if k == 0: print(base, 'sum r', r.sum(), ' 2K', 2*(np.r_[u,v] @ Kmat(X,Y) @ np.r_[u,v]), ' max r_v', r.max(), ' min', r.min())
    print(base, f'fraction of vertices with r_v > 0: {pos}/{tot}')
