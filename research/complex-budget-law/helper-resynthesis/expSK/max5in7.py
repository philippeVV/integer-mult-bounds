"""For every G26 edge L7=<p_t,p_a>^perp (26 ports): max #ports in a 5-dim subspace of L7, and max #ports in a
6-dim subspace of L7 (hyperplane). A register alone at its L6 (m=1) carries c(r,6)=c(r,5) <= max5."""
import json, numpy as np
from collections import Counter
P=json.load(open('instA.json'))['ports']; dot=lambda a,b:(a&b).bit_count()&1; n=112
g=json.load(open('g26.json')); G={int(k):v for k,v in g['G'].items()}
M=np.zeros((512,2),dtype=np.uint64)
for q in range(1,512):
    m=0
    for i,p in enumerate(P):
        if dot(p,q)==0: m|=1<<i
    M[q,0]=m&((1<<64)-1); M[q,1]=m>>64
def pc(x): return np.unpackbits(x.view(np.uint8).reshape(x.shape[0],-1),axis=1).sum(1)
c5=Counter(); c6=Counter()
for t in range(n):
    for a in G[t]:
        if a<t: continue
        base=M[P[t]]&M[P[a]]
        cs=np.array([c for c in range(1,512) if c not in (P[t],P[a],P[t]^P[a])])
        X=base&M[cs]; h=pc(X); c6[int(h.max())]+=1
        b5=0
        for c,xc in zip(cs,X):
            y=pc(xc&M[cs]); 
            # exclude d in span(p_t,p_a,c) (gives 6-dim or 7-dim, counts >= anyway); those are dim>5 so drop
            sp={0,P[t],P[a],P[t]^P[a]}; sp|={x^c for x in sp}
            mask=np.array([d not in sp for d in cs])
            v=int(y[mask].max()); b5=max(b5,v)
        c5[b5]+=1
print('edges',sum(c5.values()),' max ports in hyperplane of L7',dict(c6),' max ports in 5-dim subspace of L7',dict(c5))
