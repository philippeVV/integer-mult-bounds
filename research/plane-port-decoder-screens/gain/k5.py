import sys, os, itertools
sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','lower'))
from ports import P, dot
from pysat.solvers import Cadical153
n=len(P); A=[set(b for b in range(n) if b!=a and dot(P[a],P[b])==0) for a in range(n)]
K5=[]
for a in range(n):
    for b in A[a]:
        if b<=a: continue
        Cab=A[a]&A[b]
        for c in Cab:
            if c<=b: continue
            Cabc=Cab&A[c]
            for d in Cabc:
                if d<=c: continue
                for e in Cabc&A[d]:
                    if e>d: K5.append((a,b,c,d,e))
print('K5 count',len(K5))
for p in (2,3):
    s=Cadical153(); x=lambda t,j: t*p+j+1
    for t in range(n): s.add_clause([x(t,j) for j in range(p)])
    for k in K5:
        for j in range(p): s.add_clause([-x(t,j) for t in k])
    s.add_clause([x(0,0)])
    res=s.solve(); print('p',p,'K5-free colouring exists:',res, flush=True)
