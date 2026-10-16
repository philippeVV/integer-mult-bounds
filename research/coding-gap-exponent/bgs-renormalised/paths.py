"""Primal check restricted to hop-shortest (5-hop) paths: if rate 147/157 is routable, MCF* = 147/157 exactly
(matching the uniform dual), every edge is saturated, and the optimal face is pinned down."""
import sys; sys.path.insert(0, __file__.rsplit('/',1)[0])
from bh import CORE, SESS
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
from collections import deque
adj = [[] for _ in range(182)]
for e,(u,v) in enumerate(CORE): adj[u].append((v,e)); adj[v].append((u,e))
def bfs(s):
    d=[-1]*182; d[s]=0; dq=deque([s])
    while dq:
        x=dq.popleft()
        for y,_ in adj[x]:
            if d[y]<0: d[y]=d[x]+1; dq.append(y)
    return d
def allsp(s,t):
    ds, dt = bfs(s), bfs(t); L = ds[t]; out=[]
    def rec(x, path):
        if x==t: out.append(tuple(path)); return
        for y,e in adj[x]:
            if ds[y]==ds[x]+1 and dt[y]==L-ds[y]: rec(y, path+[e])
    rec(s,[]); return out
PATHS=[allsp(s,t) for q,s,t in SESS]
np_=[len(p) for p in PATHS]; print("5-hop paths per session: min",min(np_),"max",max(np_),"total",sum(np_))
if __name__=="__main__":
    k=len(SESS); m=len(CORE); cols=[(i,p) for i in range(k) for p in PATHS[i]]; n=len(cols)+1
    A_eq=sp.lil_matrix((k,n)); A_ub=sp.lil_matrix((m,n))
    for j,(i,p) in enumerate(cols):
        A_eq[i,j]=1
        for e in p: A_ub[e,j]+=1
    for i in range(k): A_eq[i,n-1]=-1
    c=np.zeros(n); c[-1]=-1
    res=linprog(c,A_ub=A_ub.tocsr(),b_ub=np.ones(m),A_eq=A_eq.tocsr(),b_eq=np.zeros(k),bounds=[(0,None)]*n,method='highs')
    print("restricted-path MCF rate:", -res.fun, " 147/157 =", 147/157)
    load=A_ub.tocsr()@res.x; print("edge loads: min %.9f max %.9f"%(load.min(),load.max()))
