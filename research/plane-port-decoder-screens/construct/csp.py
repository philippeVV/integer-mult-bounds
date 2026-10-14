import sys, random; sys.path.insert(0,'.')
from e8 import *
L,labs,G,fr,ports=p112()
n=len(ports); U=[labs[t] for t in ports]
odd=[[j for j in range(n) if j!=i and (ip(L[ports[i]],L[ports[j]])&1)] for i in range(n)]
dotm={}
def d(a,b): 
    k=(a,b)
    if k not in dotm: dotm[k]=gform(a,b,G)
    return dotm[k]
def okplane(u,r): return r!=0 and r!=u and (d(r,r)^d(u,r))==1
def orth(i,ri,j,rj): return d(U[i],rj)==0 and d(ri,U[j])==0 and d(ri,rj)==0
