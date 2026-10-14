import numpy as np
def r_counts(M, K):
    """Smith invariant valuations over Z/2^K (K<=40) for integer matrix M; returns list of valuations (<K)."""
    mod = 1 << K
    A = np.array(M, dtype=np.int64) % mod
    vals = []
    while A.size:
        nz = A != 0
        if not nz.any(): break
        # valuation of each entry
        low = A & -A
        v = np.where(nz, np.log2(np.where(nz, low, 1)).astype(np.int64), K)
        i, j = np.unravel_index(np.argmin(v), v.shape); vv = int(v[i, j]); vals.append(vv)
        p = int(A[i, j]) >> vv; pinv = pow(p, -1, mod)
        col = (A[:, j] >> vv) * pinv % mod          # factors
        A = (A - (col[:, None] * A[i][None, :]) % mod) % mod
        A = np.delete(np.delete(A, i, 0), j, 1)
    return vals
def r_j(M, j): return sum(1 for v in r_counts(M, j + 1) if v <= j)
