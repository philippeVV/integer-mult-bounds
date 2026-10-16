"""Singer-equivariant Tanner codes on PG(2,17): dim W = sum_chi dim(A ∩ Λ_chi B), Λ_chi = diag(ζ^{chi d})_{d∈D}.
Exact mod split primes P (dim mod P >= dim over the number field, so 'mod-P <= 1005' is a rigorous stop for that code)."""
import numpy as np, itertools, math
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from geom import singer
from modp import batch_rank, split_prime, root_of_unity, prim_root
q = 17
N, D = singer(q); n = q+1; k = n//2; E = N*n
L = math.lcm(N, q-1, q, 8, 3)
PRIMES = [split_prime(L, 10**9), split_prime(L, 2*10**9 - 10**8)]
class Ctx:
    def __init__(s, P):
        s.P = P; s.g = prim_root(P); s.z = root_of_unity(P, N, s.g)
        s.zp = [pow(s.z, e, P) for e in range(N)]
        s.Lam = np.array([[s.zp[(c*d) % N] for d in D] for c in range(N)], dtype=np.int64)  # (N, n)
def blocks(ctx, A, B):
    P = ctx.P
    A = np.asarray(A, dtype=np.int64) % P; B = np.asarray(B, dtype=np.int64) % P
    kA, kB = A.shape[0], B.shape[0]
    LB = (B[None, :, :] * ctx.Lam[:, None, :]) % P          # (N,kB,n)
    M = np.concatenate([np.broadcast_to(A, (N,)+A.shape), LB], axis=1)
    return kA + kB - batch_rank(M, P)
def rank(ctx, A): return int(batch_rank(np.asarray(A, dtype=np.int64)[None], ctx.P)[0])
_SUBS = None
def mds_bad(ctx, A, maxcheck=None):
    """number of singular k-subsets of columns (exact mod P; all C(18,9)=48620)."""
    global _SUBS
    A = np.asarray(A, dtype=np.int64) % ctx.P
    kk = A.shape[0]
    subs = np.array(list(itertools.combinations(range(A.shape[1]), kk)))
    bad = 0
    for i in range(0, len(subs), 8000):
        S = subs[i:i+8000]
        M = A[:, S].transpose(1, 0, 2)
        bad += int((batch_rank(M, ctx.P) < kk).sum())
    return bad
