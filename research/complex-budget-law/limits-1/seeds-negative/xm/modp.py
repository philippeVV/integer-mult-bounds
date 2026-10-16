"""Exact linear algebra mod a prime P < 2^31, batched over a leading axis."""
import numpy as np, math
def isprime(n): return n > 1 and all(n % r for r in range(2, int(n**.5)+1))
def split_prime(L, lo):
    P = lo - lo % L + 1
    while not isprime(P): P += L
    return P
def prim_root(P):
    fs = []; x = P-1; r = 2
    while r*r <= x:
        if x % r == 0:
            fs.append(r)
            while x % r == 0: x //= r
        r += 1
    if x > 1: fs.append(x)
    for g in range(2, P):
        if all(pow(g, (P-1)//f, P) != 1 for f in fs): return g
def root_of_unity(P, n, g=None):
    g = g or prim_root(P); assert (P-1) % n == 0
    return pow(g, (P-1)//n, P)
def batch_rank(M, P):
    """M: (b, n, m) int array; returns ranks (b,) mod P."""
    M = np.array(M, dtype=np.int64) % P
    b, n, m = M.shape
    r = np.zeros(b, dtype=np.int64); ar = np.arange(b); rows = np.arange(n)
    for c in range(m):
        mask = (rows[None, :] >= r[:, None]) & (M[:, :, c] != 0)
        has = mask.any(1)
        if not has.any(): continue
        piv = mask.argmax(1)
        idx = ar[has]; rr = r[has]; pp = piv[has]
        tmp = M[idx, rr].copy(); M[idx, rr] = M[idx, pp]; M[idx, pp] = tmp
        pv = M[idx, rr, c]
        inv = np.array([pow(int(x), P-2, P) for x in pv], dtype=np.int64)
        M[idx, rr] = (M[idx, rr] * inv[:, None]) % P
        prow = M[idx, rr]                       # (h, m)
        col = M[idx, :, c].copy()               # (h, n)
        col[np.arange(len(idx)), rr] = 0
        M[idx] = (M[idx] - (col[:, :, None] * prow[:, None, :]) % P) % P
        r[has] += 1
    return r
def batch_det_nonzero(M, P):
    b, n, _ = M.shape
    return batch_rank(M, P) == n
