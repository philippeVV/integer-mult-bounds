"""Generic-realisation screen (GRS) for Braverman-He scenario A, plus signed-label frame families.

Theorem used (telescoping of stronger-rank-screens Sec.1, without the monomial hypothesis):
for ANY linear gate assignment G over a char-0 field with the circuit's incidence (gates may be
arbitrary rectangular matrices; C-ports carry 0, J-ports are discarded -- both give zero terms),
    sum_e (b_e a_e^T) (x) (M_head(e) - M_tail(e)) = D_out (T_G (x) I) - (T_G (x) I) D_in =: Z_G,
so  s >= rank_Q Z_G,  Z_G[j,k] = T_G[j,k] (M_out(j) - M_in(k)).
Interior frames do not enter Z_G: the bound depends only on terminal frames.
We take integer G, reduce mod a prime p (p not dividing any frame denominator); then
rank_p(Z_G mod p) <= rank_Q Z_G <= s is an exact certificate.
"""
import os, sys, contextlib, io, random
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
N1 = HERE   # bh_build.py, bh_data.py (expN1), copied unchanged
sys.path.insert(0, N1)
_cwd = os.getcwd(); os.chdir(HERE)
with contextlib.redirect_stdout(io.StringIO()):
    import bh_build as bh          # re-runs the exact BH checks; writes ./bh_graph.json here
os.chdir(_cwd)

P = 2147483629                     # prime < 2^31 (products fit in int64)
B = bh.B; NB = len(B); IDX = {q: i for i, q in enumerate(B)}
PHYS = set(bh.phys)

def ports(v):
    ins = [q for q in bh.din[v] if q in PHYS or (q in bh.Rs and q not in bh.Cs and q not in bh.Js)]
    outs = [q for q in bh.dout[v] if q in PHYS or (q in bh.Rs and q not in bh.Cs and q not in bh.Js)]
    return ins, outs

PORTS = {v: ports(v) for v in range(182)}
NPAR = sum(len(i) * len(o) for i, o in PORTS.values())

def transfer(rng, p=P, gates=None):
    """T_G mod p (sinks x sources) for random integer gates (or supplied dict v->matrix)."""
    a = {}
    T = np.zeros((NB, NB), dtype=np.int64)
    for v in bh.PI:
        ins, outs = PORTS[v]
        X = np.zeros((len(ins), NB), dtype=np.int64)
        for r, q in enumerate(ins):
            if q in PHYS: X[r] = a[q]
            else: X[r, IDX[q]] = 1
        G = gates[v] if gates is not None else rng.integers(0, p, size=(len(outs), len(ins)), dtype=np.int64)
        Y = matmul(G, X, p)
        for r, q in enumerate(outs):
            if q in PHYS: a[q] = Y[r]
            else: T[IDX[q]] = Y[r]
    return T

def matmul(A, Bm, p=P):
    # exact mod-p product with int64 (split to avoid overflow in the dot accumulation)
    A = A % p; Bm = Bm % p
    lo = Bm & 0xFFFF; hi = Bm >> 16
    r1 = np.zeros((A.shape[0], Bm.shape[1]), dtype=np.int64)
    r2 = np.zeros_like(r1)
    for k in range(A.shape[1]):      # small inner dims here (<=5 or modest)
        r1 = (r1 + np.outer(A[:, k], lo[k])) % p
        r2 = (r2 + np.outer(A[:, k], hi[k])) % p
    return (r1 + (r2 * 65536) % p) % p

def rank_mod(M, p=P):
    M = M.copy() % p; n, c = M.shape; r = 0
    for col in range(c):
        if r == n: break
        nz = np.nonzero(M[r:, col])[0]
        if len(nz) == 0: continue
        piv = r + nz[0]
        if piv != r: M[[r, piv]] = M[[piv, r]]
        inv = pow(int(M[r, col]), p - 2, p)
        M[r] = (M[r] * inv) % p
        rows = np.nonzero(M[:, col])[0]; rows = rows[rows != r]
        if len(rows):
            f = M[rows, col][:, None]
            M[rows] = (M[rows] - mulmod(f, M[r][None, :], p)) % p
        r += 1
    return r

def mulmod(f, row, p=P):
    lo = row & 0xFFFF; hi = row >> 16
    return ((f * lo) % p + ((f * hi) % p) * 65536 % p) % p

def inv_mod_matrix(A, p=P):
    n = A.shape[0]; M = np.concatenate([A % p, np.eye(n, dtype=np.int64)], axis=1)
    for col in range(n):
        nz = np.nonzero(M[col:, col])[0]
        if len(nz) == 0: raise ValueError('singular mod p')
        piv = col + nz[0]; M[[col, piv]] = M[[piv, col]]
        M[col] = (M[col] * pow(int(M[col, col]), p - 2, p)) % p
        rows = np.nonzero(M[:, col])[0]; rows = rows[rows != col]
        if len(rows): M[rows] = (M[rows] - mulmod(M[rows, col][:, None], M[col][None, :], p)) % p
    return M[:, n:]

def frac_mod(num, den, p=P):
    return (num % p) * pow(den % p, p - 2, p) % p

def Zmat(T, Min, Mout, p=P):
    """Z[j,k] = T[j,k] (Mout[j] - Min[k]); Min/Mout: lists of m x m int64 mod p."""
    W = T.shape[0]; m = Min[0].shape[0]
    Z = np.zeros((W * m, W * m), dtype=np.int64)
    for j in range(W):
        for k in np.nonzero(T[j])[0]:
            Z[j*m:(j+1)*m, k*m:(k+1)*m] = mulmod(np.int64(T[j, k]), (Mout[j] - Min[k]) % p, p)
    return Z
