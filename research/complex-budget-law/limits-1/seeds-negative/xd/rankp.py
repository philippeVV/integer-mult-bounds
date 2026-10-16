"""Blocked rank mod a prime p with float32 BLAS updates (entries < p, panel width w: w*p^2 < 2^24)."""
import numpy as np, numba as nb
@nb.njit(cache=True)
def _panel(P, p, inv, piv):
    # P: m x w int64 copy; returns rank r, pivot row ids in piv[:r], RREF rows stored in out (r x w), pivot cols
    m, w = P.shape; r = 0
    rows = np.arange(m); pc = np.zeros(w, dtype=np.int64)
    for c in range(w):
        q = -1
        for i in range(r, m):
            if P[i, c] % p != 0: q = i; break
        if q < 0: continue
        if q != r:
            for k in range(w):
                t = P[r, k]; P[r, k] = P[q, k]; P[q, k] = t
            t = rows[r]; rows[r] = rows[q]; rows[q] = t
        iv = inv[P[r, c] % p]
        for k in range(w): P[r, k] = (P[r, k] * iv) % p
        for i in range(m):
            if i != r and P[i, c] % p != 0:
                f = P[i, c] % p
                for k in range(w): P[i, k] = (P[i, k] - f * P[r, k]) % p
        pc[r] = c; r += 1
    for i in range(r): piv[i] = rows[i]
    return r, pc
def rank_modp(A, p, w=None):
    A = np.asarray(A, dtype=np.float32) % p
    if w is None: w = max(8, min(512, (1 << 24) // (p * p) // 2))
    inv = np.zeros(p, dtype=np.int64)
    for a in range(1, p): inv[a] = pow(a, p - 2, p)
    rank = 0
    while A.shape[0] and A.shape[1]:
        ww = min(w, A.shape[1]); P = A[:, :ww].astype(np.int64)
        piv = np.zeros(ww, dtype=np.int64)
        P0 = P.copy(); r, pc = _panel(P, p, inv, piv)
        if r == 0: A = A[:, ww:]; continue
        R = piv[:r]
        # solve: every row i of P0 = x_i @ B where B = rows of RREF basis (P[:r]); x_i = P0[i, pc[:r]]
        Bpan = P[:r]  # RREF of P0[R]; T with T @ P0[R] = Bpan, T = inverse of P0[R][:,pc]
        Mr = P0[R][:, pc[:r]] % p  # r x r invertible
        Tm = inv_modp(Mr, p, inv)   # T = Mr^{-1}
        X = (P0[:, pc[:r]] % p).astype(np.float32)  # coefficient of each row on basis T@A[R]
        keep = np.ones(A.shape[0], bool); keep[R] = False
        trail = A[:, ww:]
        if trail.shape[1]:
            TB = mm(Tm.astype(np.float32), trail[R], p)       # r x n'
            upd = mm(X[keep], TB, p)
            A = (trail[keep] - upd) % p
        else:
            A = trail[keep]
        rank += r
    return rank
def mm(X, Y, p, blk=None):
    # exact product mod p in float32, chunking the inner dimension
    k = X.shape[1]; step = max(1, (1 << 24) // ((p - 1) ** 2) - 1)
    out = np.zeros((X.shape[0], Y.shape[1]), dtype=np.float32)
    for s in range(0, k, step):
        out += X[:, s:s + step] @ Y[s:s + step]
        out %= p
    return out
def inv_modp(M, p, inv):
    n = M.shape[0]; A = np.concatenate([M % p, np.eye(n, dtype=np.int64)], axis=1)
    piv = np.zeros(2 * n, dtype=np.int64); r, pc = _panel(A, p, inv, piv)
    assert r == n
    return A[:, n:]
