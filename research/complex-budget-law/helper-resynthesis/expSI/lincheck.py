"""R1 check (exact, mod p=2^61-1): with the edge multiplicities of coding.py, assign explicit linear forms in the 112
source symbols to every edge and verify (i) every edge form lies in the span of the forms entering its tail node,
(ii) every target top / centre frame can form its exact demand (alpha_t, centre A-content) from its arrivals.
Uses random substitution: forms are evaluated as exact vectors mod p; spans by Gaussian elimination mod p."""
import json, random
from fractions import Fraction as F
from collections import defaultdict
p=(1<<61)-1; md=lambda q: (q.numerator*pow(q.denominator,-1,p))%p
d=json.load(open('quads.json')); V=d['V']; C=d['C']; rects=[(tuple(A),tuple(B)) for A,B in d['rects']]
al=json.load(open('alpha.json')); I=json.load(open('instA.json')); n=112
src=open('design.py').read(); exec(src[src.index('import json'):src.index('seed=int')]); exec(src[src.index('def perp'):src.index('# ---- nodes')])
def coef(b,s):
    if b<112: return F(al[str(b)].get(str(s),'0'))
    return F(I['centres'][b-112]['src'].get(str(s),0))
def reduce(B,v):
    v=list(v)
    for piv,row in B:
        if v[piv]: c=v[piv]; v=[(x-c*y)%p for x,y in zip(v,row)]
    return v
def addb(B,v):
    v=reduce(B,v); i=next((k for k,x in enumerate(v) if x),None)
    if i is None: return False
    inv=pow(v[i],-1,p); v=[x*inv%p for x in v]
    B[:]=[(q,[(x-r[i]*y)%p for x,y in zip(r,v)]) for q,r in B]; B.append((i,v)); return True
inspan=lambda B,v: not any(reduce(B,v))
def vec(dct): 
    v=[0]*n
    for s,c in dct.items(): v[s]=(v[s]+md(F(c)))%p
    return v
Mab={(A,B):[[coef(b,s) for s in A] for b in B] for A,B in rects}
ga=defaultdict(list); sc=defaultdict(list)
for A,B in rects: ga[A].append(B); sc[B].append(A)
ok=True; edges=0
# gather: pair node P emits a basis of proj_P(R_A) to quad A; quad A emits a basis of rowspace(M_AB) to rectangle (A,B)
for A,Bs in ga.items():
    qin=[]
    for pr in qpair[A]:
        pin=[]; [addb(pin,vec({a:1})) for a in pr]                     # line registers bring x_a, x_b
        proj=[]
        for B in Bs:
            for row in Mab[(A,B)]:
                v=vec({s:c for s,c in zip(A,row) if s in pr})
                if addb(proj,v): ok&=inspan(pin,v); edges+=1; addb(qin,v)
    for B in Bs:
        rs=[]
        for row in Mab[(A,B)]:
            v=vec(dict(zip(A,row)))
            if addb(rs,v): ok&=inspan(qin,v); edges+=1
print('gather side: every emitted form lies in its node in-span:',ok)
# scatter: dual quad B in-span = rowspaces from its rectangles; emits to D a basis of {row_b: b in D}; D emits row_D(b) to top b
arr=defaultdict(list); ok2=True
for B,As in sc.items():
    bin_=[]
    for A in As:
        for row in Mab[(A,B)]: addb(bin_,vec(dict(zip(A,row))))
    for pr in qpair[B]:
        for b in pr:
            v=vec({s:c for A in As for s,c in zip(A,Mab[(A,B)][B.index(b)])})
            ok2&=inspan(bin_,v); arr[(b,pr)].append(v)
# D -> top b carries the sum over its dual quads; top sums its arrivals
good=0
for b in range(120):
    tot=[0]*n
    for (bb,pr),vs in arr.items():
        if bb==b:
            for v in vs: tot=[(x+y)%p for x,y in zip(tot,v)]
    want=vec({s:coef(b,s) for s in range(n) if coef(b,s)!=0}) if b<112 else vec({int(k):c for k,c in I['centres'][b-112]['src'].items()})
    if b>=112: want=[0 if not (dot(V[s],V[b])==0) else w for s,w in enumerate(want)] if False else want
    good+= tot==want
print('scatter side forms in span:',ok2,' tops/centres whose arrivals sum exactly to their demand:',good,'/ 120')
