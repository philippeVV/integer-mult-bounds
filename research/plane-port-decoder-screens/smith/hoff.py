import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla, sys
from omega import *
D = np.array([[dot(a, b) for b in range(512)] for a in range(512)], dtype=np.uint8)
Q = np.array([p[0] for p in planes]); R = np.array([p[1] for p in planes])
n = len(planes)
rows = []; cols = []
rk2rows = []; rk2cols = []
for i in range(n):
    a = D[Q[i], Q]; b = D[Q[i], R]; c = D[R[i], Q]; d = D[R[i], R]
    orth = np.nonzero((a | b | c | d) == 0)[0]
    rows.append(np.full(len(orth), i)); cols.append(orth)
    r1 = np.nonzero(((a & d) ^ (b & c)) == 0)[0]   # det == 0 -> rank <= 1 (incl. 0)
    r1 = r1[r1 != i]
    rk2rows.append(np.full(len(r1), i)); rk2cols.append(r1)
for name, rr, cc in (('orth', rows, cols), ('rank<=1', rk2rows, rk2cols)):
    A = sp.csr_matrix((np.ones(sum(len(c) for c in cc), dtype=np.float64), (np.concatenate(rr), np.concatenate(cc))), shape=(n, n))
    deg = int(A[0].sum())
    lmin = sla.eigsh(A, k=1, which='SA', return_eigenvectors=False)[0]
    lmax = sla.eigsh(A, k=1, which='LA', return_eigenvectors=False)[0]
    hb = n * (-lmin) / (deg - lmin)
    print(name, 'deg', deg, 'lmax', round(lmax, 4), 'lmin', round(lmin, 4), 'Hoffman alpha <=', round(hb, 2), 'density <=', round(hb / n, 5), '-> x112 =', round(112 * hb / n, 2), flush=True)
    sp.save_npz(f'omega_{name.replace("<=","le")}.npz', A)
