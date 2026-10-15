import sys as _s, os as _o; _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "..", "helper-bounds"))  # HB_PATH
"""Build the quad hub-tree design on the #1094 interface (relaxed OR-reach model), decompose into register pieces,
assign source/target/helper roles (max-weight matching), write design.json. Verification is in verify.py."""
import json, sys, random, math, itertools
from collections import Counter, defaultdict
sys.path.insert(0,'../expBV'); from wordmodel import basis
import numpy as np
from scipy.optimize import linear_sum_assignment
seed=int(sys.argv[1]) if len(sys.argv)>1 else 0; random.seed(seed)
d=json.load(open('quads.json')); V=d['V']; C=d['C']; rects=[(tuple(A),tuple(B)) for A,B in d['rects']]
n=112; dot=lambda a,b:(a&b).bit_count()&1
def perp(B): return basis([x for x in range(1,512) if all(dot(x,b)==0 for b in B)])
def span(B):
    s=[0]
    for v in B: s=s+[x^v for x in s]
    return s
def contained(a,b): return basis(list(a)+list(b))==tuple(b)
cls={x:i for i,c in enumerate(C) for x in c}
# ---- per-class pairing choice (Fano: e1 on 3 lines, e2,e3 on 2 lines each, e1,e2,e3 collinear)
def pairing_choice(c):
    X=basis([V[x] for x in c]); o=V[c[0]]
    E=[x for x in span(X) if x and not x.bit_count()&1]          # 7 even nonzero
    lines=[frozenset(L) for L in {frozenset((a,b,a^b)) for a in E for b in E if a!=b}]
    best=None
    for e1 in E:
        for L0 in [L for L in lines if e1 in L]:
            e2,e3=sorted(L0-{e1})
            ch={}
            for L in lines:
                if e1 in L: ch[L]=e1
                elif e2 in L: ch[L]=e2
                else: ch[L]=e3
            if all(ch[L] in L for L in lines): best=ch; break
        if best: break
    return best
def quadkey(Q): return basis([V[x] for x in Q])
qpair={}    # quad (tuple of idx) -> list of 2 pairs (tuples of idx)
for c in C:
    ch=pairing_choice(c); vid={V[x]:x for x in c}
    E2s={}
    for L,e in ch.items():
        for x in c:
            Q=tuple(sorted(vid[V[x]^y] for y in [0,*L]))
            if Q not in qpair:
                prs=sorted({tuple(sorted((a,vid[V[a]^e]))) for a in Q}); assert len(prs)==2
                qpair[Q]=prs
# ---- nodes (frames) and edge multiset
G=defaultdict(int)   # (u,v) -> multiplicity ; nodes are basis tuples
role_line={}; role_top={}
L=lambda a:(V[a],)
PN=lambda pr: basis([V[pr[0]],V[pr[1]]])
TOP=lambda b: perp([V[b]])
DP=lambda pr: perp([V[pr[0]],V[pr[1]]])
ga=Counter(A for A,B in rects); sc=Counter(B for A,B in rects)
for A in ga:
    for pr in qpair[A]:
        G[(PN(pr),quadkey(A))]+=1
        for a in pr: G[(L(a),PN(pr))]=1          # one line register per used pair node
for B in sc:
    for pr in qpair[B]:
        G[(perp(quadkey(B)),DP(pr))]+=1
        for b in pr: G[(DP(pr),TOP(b))]=1
for i,(A,B) in enumerate(rects):
    QA=quadkey(A); QBp=perp(quadkey(B)); sp=set(span(QA))
    cand=[x for x in span(QBp) if x not in sp]; x=random.choice(cand); X=basis(list(QA)+[x])
    sx=set(span(X)); y=random.choice([z for z in span(QBp) if z not in sx]); Y=basis(list(X)+[y])
    assert contained(QA,X) and contained(X,Y) and contained(Y,QBp)
    G[(QA,X)]+=1; G[(X,Y)]+=1; G[(Y,QBp)]+=1
# ---- piece decomposition: at each node pair in-edges with out-edges at random
inn=defaultdict(list); out=defaultdict(list)
eid=[]
for (u,v),m in G.items():
    for _ in range(m): k=len(eid); eid.append((u,v)); out[u].append(k); inn[v].append(k)
nxt={}
nodes=set(inn)|set(out)
for N in nodes:
    a=inn[N][:]; b=out[N][:]; random.shuffle(a); random.shuffle(b)
    for i,j in zip(a,b): nxt[i]=j
hasprev=set(nxt.values())
pieces=[]
for k in range(len(eid)):
    if k in hasprev: continue
    ch=[eid[k][0],eid[k][1]]; j=k
    while j in nxt: j=nxt[j]; ch.append(eid[j][1])
    pieces.append(ch)
dimc=Counter(len(N) for p in pieces for N in p)
print('pieces',len(pieces),'registers per dim',dict(sorted(dimc.items())))
# ---- roles: sources (piece starting at own line), targets (piece ending at own top); max weight matching
f=lambda w: w*math.log(73/w)
lineof={L(a):a for a in range(n)}; topof={TOP(b):b for b in range(n)}
slots=[('S',a) for a in range(n)]+[('T',b) for b in range(n)]
W=np.zeros((len(pieces),len(slots)))
for i,p in enumerate(pieces):
    w=f(len(p[-1])-len(p[0]))
    if p[0] in lineof: W[i,lineof[p[0]]]=w
    if p[-1] in topof: W[i,n+topof[p[-1]]]=w
r,c=linear_sum_assignment(-W)
role=['H']*len(pieces)
for i,j in zip(r,c):
    if W[i,j]>0: role[i]=slots[j]
assert sum(1 for x in role if x!='H')==224, 'not all data roles filled'
hw=Counter(len(p[-1])-len(p[0]) for p,x in zip(pieces,role) if x=='H')
mass=sum(k*f(w) for w,k in hw.items())
print('helpers',sum(hw.values()),'width inventory',dict(sorted(hw.items())),'log-mass %.1f'%mass)
json.dump(dict(pieces=[[list(N) for N in p] for p in pieces],role=role,inventory=hw,mass=mass,seed=seed),open(f'design{seed}.json','w'))
