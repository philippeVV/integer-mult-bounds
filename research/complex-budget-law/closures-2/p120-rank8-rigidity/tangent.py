"""Exact tangent-space certificate for rank-8 rigidity at N0 = G/2 on the E8 P120 labels.
Variety V = {(X,Y) in Q^{120x8} x Q^{120x8} : X_a.Y_b = delta_ab for all (a,b) with a=b or a _|_ b}.
At (X0,Y0) = (r/1, r/2) [r = roots, norm 2; X0 Y0^T = G/2], compute rank of the Jacobian J (7680 x 1920) exactly
(integer matrix, FLINT fmpz_mat rank and nmod_mat rank) and the rank of the 183 gauge vectors
(GL8: dX = X0 M, dY = -Y0 M^T; diagonal: dX = D X0, dY = -D Y0) and check J*gauge = 0.
If rank J = 1920 - 183, the Zariski tangent space equals the tangent space of the gauge orbit, so V is smooth of
dimension 183 at (X0,Y0) and coincides near it with the orbit {(D X0 M, D^-1 Y0 M^-T)}: G/2 is locally rigid."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, flint, time
from e8 import R, G, n
K=8
X0=R.astype(np.int64)        # doubled coords: X0 Y0^T = R R^T / 8 *? -> use integer scaling, same tangent space
Y0=R.astype(np.int64)        # X0 Y0^T = 4G (scalar multiple of G/2 -> same variety up to scaling of equations)
cons=[(a,b) for a in range(n) for b in range(n) if a==b or G[a,b]==0]
m=len(cons); nv=2*n*K
J=np.zeros((m,nv),dtype=np.int64)
for e,(a,b) in enumerate(cons):
    J[e,a*K:(a+1)*K]+=Y0[b]          # d/dX_a
    J[e,n*K+b*K:n*K+(b+1)*K]+=X0[a]  # d/dY_b
gauge=[]
for i in range(K):
    for j in range(K):
        M=np.zeros((K,K),dtype=np.int64); M[i,j]=1
        gauge.append(np.concatenate([(X0@M).ravel(),(-(Y0@M.T)).ravel()]))
for a in range(n):
    D=np.zeros((n,1),dtype=np.int64); D[a]=1
    gauge.append(np.concatenate([(D*X0).ravel(),(-(D*Y0)).ravel()]))
Gm=np.array(gauge)
print('equations',m,'unknowns',nv)
assert not (J@Gm.T).any(); print('J * gauge = 0 : OK')
t=time.time()
rg=flint.fmpz_mat(Gm.tolist()).rank(); print('gauge rank (exact)',rg)
p=2**61-1
rJp=flint.nmod_mat(J.tolist(),p).rank(); print('rank J mod p =',rJp,'(%.0fs)'%(time.time()-t))
# rank over Q >= rank mod p, and rank over Q <= nv - rg (gauge in kernel): equality certifies
print('nv - gauge rank =',nv-rg)
if rJp==nv-rg: print('PASS_TANGENT: rank_Q(J) = %d exactly; tangent space = gauge orbit tangent (dim %d)'%(rJp,rg))
else: print('FAIL: gap',nv-rg-rJp)
