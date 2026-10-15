import sys as _s, os as _o; _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "..", "helper-bounds"))  # HB_PATH
import json,sys
from collections import Counter,defaultdict
sys.path.insert(0,'../expBV'); from wordmodel import basis
P=json.load(open('instA.json'))['ports']; cen=json.load(open('cen.json'))['cen']; n=112
dot=lambda a,b:(a&b).bit_count()&1; isp={p:i for i,p in enumerate(P)}
cl=defaultdict(list)
for s in range(n): cl[tuple(cen[s])].append(s)
print('classes',len(cl),Counter(len(v) for v in cl.values()))
T=lambda Q: frozenset(i for i,p in enumerate(P) if all(dot(p,q)==0 for q in Q))
def span(B):
    s=[0]
    for v in B: s=s+[x^v for x in s]
    return s
C=list(cl.values())
for c in C[:3]:
    B=basis([P[i] for i in c]); sp=span(B)
    print('class',c,'dim',len(B),'ports in span',sum(x in isp for x in sp),'odd in span',sum(bin(x).count('1')%2 for x in sp),
      'internal dots',Counter(dot(P[a],P[b]) for a in c for b in c if a<b),'pair-level sizes',Counter(len(T([P[a],P[b]])) for a in c for b in c if a<b))
# orthogonality between classes: for class pairs, count of orthogonal port pairs
M=Counter()
for i,a in enumerate(C):
    for j,b in enumerate(C):
        M[(i,j)]=sum(dot(P[x],P[y])==0 for x in a for y in b)
print('class-pair orth counts',Counter(M.values()))
print('row of class0',[M[(0,j)] for j in range(14)])
