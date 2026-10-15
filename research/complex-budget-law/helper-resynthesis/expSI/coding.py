import sys as _s, os as _o; _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "..", "helper-bounds"))  # HB_PATH
"""Coding-aware version of the quad design (first lifting refinements):
 (R1) linear capacity: every register edge carries ONE linear form; a node can emit only forms in the span of what entered.
      Edge multiplicities are set from the exact coefficient matrices (alpha for target tops, centre A-content for U_q).
 (R2) audit hygiene 'non-target reads target': a target can only absorb, so a target register may carry only its final
      edge D -> p_t^perp (it reads the registers at D); every relay edge is a source or a helper.
Writes coding{seed}.json in design format (pieces/role) for verify.py."""
import json, sys, random, math
from fractions import Fraction as F
from collections import Counter, defaultdict
sys.path.insert(0,'../expBV'); from wordmodel import basis
seed=int(sys.argv[1]) if len(sys.argv)>1 else 0; TROLE=sys.argv[2] if len(sys.argv)>2 else 'final'
random.seed(seed)
exec(open('design.py').read().split('# ---- nodes')[0].split('seed=int')[0])   # imports only
d=json.load(open('quads.json')); V=d['V']; C=d['C']; rects=[(tuple(A),tuple(B)) for A,B in d['rects']]
n=112; dot=lambda a,b:(a&b).bit_count()&1
src=open('design.py').read(); exec(src[src.index('def perp'):src.index('# ---- nodes')])
al=json.load(open('alpha.json')); I=json.load(open('instA.json'))
def coef(b,s):
    if b<112: return F(al[str(b)].get(str(s),'0'))
    return F(I['centres'][b-112]['src'].get(str(s),0))
def rank(M):
    M=[list(r) for r in M if any(r)]; r=0
    if not M: return 0
    for c in range(len(M[0])):
        p=next((i for i in range(r,len(M)) if M[i][c]!=0),None)
        if p is None: continue
        M[r],M[p]=M[p],M[r]
        for i in range(len(M)):
            if i!=r and M[i][c]!=0:
                f=M[i][c]/M[r][c]; M[i]=[x-f*y for x,y in zip(M[i],M[r])]
        r+=1
    return r
Mab={(A,B):[[coef(b,s) for s in A] for b in B] for A,B in rects}
rk={r:rank(m) for r,m in Mab.items()}
print('rectangle ranks', Counter((B[0]>=112, k) for (A,B),k in rk.items()))
G=defaultdict(int)
L=lambda a:(V[a],); PN=lambda pr: basis([V[pr[0]],V[pr[1]]]); TOP=lambda b: perp([V[b]]); DP=lambda pr: perp([V[pr[0]],V[pr[1]]])
ga=defaultdict(list); sc=defaultdict(list)
for A,B in rects: ga[A].append(B); sc[B].append(A)
for A,Bs in ga.items():
    rows=[r for B in Bs for r in Mab[(A,B)]]
    for pr in qpair[A]:
        idx=[A.index(x) for x in pr]
        k=rank([[r[i] for i in idx] for r in rows])        # forms needed from this pair node
        G[(PN(pr),quadkey(A))]+=k
        for a in pr: G[(L(a),PN(pr))]=1
for B,As in sc.items():
    for pr in qpair[B]:
        rows=[[x for A in As for x in Mab[(A,B)][B.index(b)]] for b in pr]
        G[(perp(quadkey(B)),DP(pr))]+=rank(rows)
        for b in pr: G[(DP(pr),TOP(b))]=1
for (A,B),k in rk.items():
    QA=quadkey(A); QBp=perp(quadkey(B)); sp_=set(span(QA))
    x=random.choice([x for x in span(QBp) if x not in sp_]); X=basis(list(QA)+[x])
    sx=set(span(X)); y=random.choice([z for z in span(QBp) if z not in sx]); Y=basis(list(X)+[y])
    G[(QA,X)]+=k; G[(X,Y)]+=k; G[(Y,QBp)]+=k
# pair nodes need >= forms they emit: in-span of a pair node is <= #line registers entering (here 2) -> fine (k<=2).
# ---- target edges: one D->top_t edge per target (R2) ; random choice among its 3
topof={TOP(b):b for b in range(n)}; lineof={L(a):a for a in range(n)}
eid=[]; inn=defaultdict(list); out=defaultdict(list)
for (u,v),m in G.items():
    for _ in range(m): k=len(eid); eid.append((u,v)); out[u].append(k); inn[v].append(k)
tedge={}
for t in range(n):
    cand=inn[TOP(t)]; tedge[t]=random.choice(cand)
tset=set(tedge.values())
nxt={}
for N in set(inn)|set(out):
    a=[k for k in inn[N] if k not in tset]; b=[k for k in out[N] if k not in tset]; random.shuffle(a); random.shuffle(b)
    for i,j in zip(a,b): nxt[i]=j
hasprev=set(nxt.values()); pieces=[]; role=[]
for k in range(len(eid)):
    if k in tset:
        pieces.append([eid[k][0],eid[k][1]]); role.append(('T',topof[eid[k][1]])); continue
    if k in hasprev: continue
    ch=[eid[k][0],eid[k][1]]; j=k
    while j in nxt: j=nxt[j]; ch.append(eid[j][1])
    pieces.append(ch); role.append('H')
for a in range(n):   # source = longest piece starting at its line
    c=[i for i,p in enumerate(pieces) if role[i]=='H' and p[0]==L(a)]
    i=max(c,key=lambda i:len(pieces[i][-1])); role[i]=('S',a)
f=lambda w: w*math.log(73/w)
hw=Counter(len(p[-1])-len(p[0]) for p,r in zip(pieces,role) if r=='H')
mass=sum(k*f(w) for w,k in hw.items())
print('edges per dim-step',dict(sorted(Counter(len(u) for u,v in eid).items())))
print('seed',seed,'helpers',sum(hw.values()),'inventory',dict(sorted(hw.items())),'log-mass %.1f'%mass)
json.dump(dict(pieces=[[list(N) for N in p] for p in pieces],role=role,inventory=hw,mass=mass,seed=seed),open(f'coding{seed}.json','w'))
