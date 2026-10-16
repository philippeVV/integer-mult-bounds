"""Rational lower-bound certificate for the core MCF (rate ~147/157) from the restricted 5-hop-path LP."""
import sys; sys.path.insert(0, __file__.rsplit('/',1)[0])
import numpy as np, scipy.sparse as sp, json
from scipy.optimize import linprog
from paths import PATHS
from bh import CORE, SESS
k=len(SESS); m=len(CORE); cols=[(i,p) for i in range(k) for p in PATHS[i]]; n=len(cols)+1
A_eq=sp.lil_matrix((k,n)); A_ub=sp.lil_matrix((m,n))
for j,(i,p) in enumerate(cols):
    A_eq[i,j]=1
    for e in p: A_ub[e,j]+=1
for i in range(k): A_eq[i,n-1]=-1
c=np.zeros(n); c[-1]=-1
res=linprog(c,A_ub=A_ub.tocsr(),b_ub=np.ones(m),A_eq=A_eq.tocsr(),b_eq=np.zeros(k),bounds=[(0,None)]*n,method='highs')
out=[(i,list(p),float(res.x[j])) for j,(i,p) in enumerate(cols) if res.x[j]>0]
json.dump(out,open('core_flow.json','w')); print("saved",len(out),"positive path flows; LP rate",-res.fun)
