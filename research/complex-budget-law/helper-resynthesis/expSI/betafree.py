import sys as _s, os as _o; _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "..", "helper-bounds"))  # HB_PATH
"""L1/L2 with beta free: alpha_t(beta) = f_t - sum_t' beta_tt' c_t' (x-parts; f_t = final x-part of target t,
c_t' = post-cut x-part of target t'). Question: can some beta put alpha_t's support inside a proper subspace of
p_t^perp (which would remove the 'second top arrival' requirement for t)? Exact test mod a large prime, for each t and
each 7-dim L7 < p_t^perp: solve  (f_t - C beta)|_{sources outside L7} = 0."""
import sys, json; sys.path.insert(0,'../expBV')
from wordmodel import Word
from fractions import Fraction as F
from collections import Counter
W=Word('../../helper-bounds/sel1094.json.gz'); w=W.w; v=W.v; n=W.n; P=w['ports']; p=2**31-1
dot=lambda a,b:(a&b).bit_count()&1
rows=[{k:F(1)} for k in range(n)]
def ap(e):
    if e[0]=='add':
        t,s=e[1],e[2]; q=F(e[3],e[4])
        for k,x in list(rows[s].items()):
            a=rows[t].get(k,0)+q*x
            if a: rows[t][k]=a
            else: rows[t].pop(k,None)
    else:
        q=F(e[2],e[3]); rows[e[1]]={k:q*x for k,x in rows[e[1]].items()}
for sec in ('precondition','compensation','A','cut'): 
    for e in w[sec]: ap(e)
C=[{k:x for k,x in rows[v+t].items() if k<v} for t in range(v)]
Cy=[{k:x for k,x in rows[v+t].items() if v<=k<2*v} for t in range(v)]
for sec in ('B','restore'):
    for e in w[sec]: ap(e)
Fx=[{k:x for k,x in rows[v+t].items() if k<v} for t in range(v)]
print('post-cut target x-support sizes',Counter(len(c) for c in C), ' y-part sizes',Counter(len(c) for c in Cy))
print('final x-support sizes',Counter(len(f) for f in Fx))
md=lambda x: (x.numerator*pow(x.denominator,-1,p))%p
import numpy as np
Cm=np.array([[md(C[t].get(s,F(0))) for t in range(v)] for s in range(v)],dtype=object)  # rows sources, cols targets
def solvable(A,b):
    # Gaussian elimination mod p on augmented matrix
    M=[list(r)+[bb] for r,bb in zip(A,b)]; r=0; ncol=len(M[0])-1
    for c in range(ncol):
        piv=next((i for i in range(r,len(M)) if M[i][c]%p),None)
        if piv is None: continue
        M[r],M[piv]=M[piv],M[r]; inv=pow(M[r][c],-1,p)
        M[r]=[x*inv%p for x in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]%p:
                f=M[i][c]; M[i]=[(x-f*y)%p for x,y in zip(M[i],M[r])]
        r+=1
    return all(any(x%p for x in row[:-1]) or row[-1]%p==0 for row in M)
res=Counter(); 
for t in range(v):
    f=[md(Fx[t].get(s,F(0))) for s in range(v)]
    perp=[x for x in range(1,512) if dot(x,P[t])==0]
    hyps=set()
    for c in range(1,512):
        if c==P[t]: continue
        hyps.add(frozenset(s for s in range(v) if dot(P[s],P[t])==0 and dot(P[s],c)==0))
    ok=0
    for H in hyps:
        out=[s for s in range(v) if s not in H]
        if solvable([list(Cm[s]) for s in out],[f[s] for s in out]): ok+=1; break
    res[ok]+=1
    if t%28==0: print('t',t,'some beta confines alpha_t to a 7-dim L7:',bool(ok),flush=True)
print('targets for which SOME beta confines alpha_t to a proper subspace:',dict(res))
