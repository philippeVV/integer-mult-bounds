"""Coefficient structure of alpha on the 4x4 rectangles: rank of alpha|AxB and, per gather quad, the dim of the span
of the s-combinations it must emit (one per rectangle), vs the 2 registers a quad node receives."""
import json
from fractions import Fraction as F
from collections import Counter, defaultdict
def rank(M):
    M=[list(r) for r in M]; r=0
    for c in range(len(M[0]) if M else 0):
        p=next((i for i in range(r,len(M)) if M[i][c]!=0),None)
        if p is None: continue
        M[r],M[p]=M[p],M[r]
        for i in range(len(M)):
            if i!=r and M[i][c]!=0:
                f=M[i][c]/M[r][c]; M[i]=[x-f*y for x,y in zip(M[i],M[r])]
        r+=1
    return r
d=json.load(open('quads.json')); rects=[(tuple(A),tuple(B)) for A,B in d['rects']]
al=json.load(open('alpha.json')); a=lambda t,s: F(al[str(t)].get(str(s),'0'))
rk=Counter(); need=defaultdict(list); vals=Counter()
for A,B in rects:
    if B[0]>=112: continue
    M=[[a(t,s) for s in A] for t in B]; r=rank(M); rk[r]+=1
    for row in M: vals[str(row)]+=0
    need[A]+= M          # quad A must make every row (t's s-combination) available downstream
print('rank of alpha restricted to port rectangles:',dict(rk))
print('gather quad: (#rows needed, dim of their span):',dict(Counter((len(v),rank(v)) for v in need.values())))
print('alpha values on rectangles',Counter(str(a(t,s)) for A,B in rects if B[0]<112 for s in A for t in B))
