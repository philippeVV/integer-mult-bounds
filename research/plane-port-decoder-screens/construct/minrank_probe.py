"""FLOAT PROBE ONLY: alternating projections for min rank completion of plane-port pattern (diag=1, required=0, free=any)."""
import sys, json; sys.path.insert(0,'.')
from csp import *
import numpy as np
R=json.load(open(sys.argv[1]))['R']; K=int(sys.argv[2])
free=np.zeros((n,n),bool)
for i in range(n):
    for j in odd[i]:
        if orth(i,R[i],j,R[j]): free[i,j]=True
rng=np.random.default_rng(0)
M=rng.standard_normal((n,n))
for it in range(3000):
    M[~free]=0; np.fill_diagonal(M,1)
    u,s,vt=np.linalg.svd(M); X=(u[:,:K]*s[:K])@vt[:K]
    err=np.linalg.norm((X-M)[~free]) 
    M=X.copy()
print('K',K,'free entries',free.sum(),'residual',err)
